import os
import numpy as np
import matplotlib.pyplot as plt
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import CoolProp.CoolProp as CP

# Initialize REFPROP/CoolProp
try:
    RP = REFPROPFunctionLibrary('C:\\Program Files\\REFPROP') # Modify your install location if necessary
    RP.SETPATHdll(os.environ.get('RPPREFIX', r"C:\Program Files\REFPROP")) # Modify your install location if necessary
    REFPROP = True
    print("Using REFPROP.")
except ValueError:
    REFPROP = False
    print("Using CoolProp.")

def PropsSI_auto(output: str, key1: str, val1: float, key2: str, val2: float, fluid: str):
    """
    Selects a fluid EOS solver depending if you have a REFPROP license. Otherwise, CoolProp will be used.
    """
    if REFPROP:
        if output == "Q":
            result = RP.REFPROPdll(
            fluid,
            key1 + key2,
            "QMASS",
            RP.MASS_BASE_SI,  # base SI units
            0, 0,             # iFlag, iUnits
            val1, val2,
            [1.0]              # composition (pure fluid)
            )
        else:
            result = RP.REFPROPdll(
            fluid,
            key1 + key2,
            output,
            RP.MASS_BASE_SI,  # base SI units
            0, 0,             # iFlag, iUnits
            val1, val2,
            [1.0]              # composition (pure fluid)
            )

        return result.Output[0]
    else:
        return  CP.PropsSI(output, key1, val1, key2, val2, fluid)


class Node():
    """
    Node class. State defined by total density d (kg/m^3) and enthalpy K(J).
    Initialized by fluid, mass m (kg), volume V (L), tempurature T (K), and name.
    """
    def __init__(self, fluid, m, V, T, name="node"):
        self.fluid = fluid
        self.m = float(m) # node mass [kg]
        self.V = float(V) / 1000.0  # convert L -> m^3
        self.name = name

        # Initialize state using T and density computed from m/V
        self.d = self.m / self.V
        # specific enthalpy from (D,T)
        h_spec = PropsSI_auto('H', 'D', self.d, 'T', float(T), self.fluid)  # J/kg
        self.H = self.m * h_spec # total enthalpy [J]
        # derived (will also populate m_l, m_v)
        self._flash_from_DH(self.d, self.H)

        # node history dict initialization
        self.history = {k: [] for k in ["time","Q","P","T","H","h","d","m","m_l","m_v", "fill_level"]}

    def _flash_from_DH(self, d, H):
        """
        Given bulk density d (kg/m3) and total enthalpy H (J),
        compute T, P, h and split m into m_l, m_v if two-phase.
        """
        m = self.m
        if m <= 0:
            # safety floor
            m = 1e-12
            self.m = m

        h = H / m  # specific enthalpy J/kg
        # try to get T and P from (D,H)
        try:
            T = PropsSI_auto('T', 'D', d, 'H', h, self.fluid)
            P = PropsSI_auto('P', 'D', d, 'H', h, self.fluid)
            phase = CP.PhaseSI('D', d, 'H', h, self.fluid) # use only CoolProp here, REFPROP phase lookup behaves weirdly
        except Exception as e:
            raise RuntimeError(f"CoolProp lookup failed in flash: d={d}, h={h}, err={e}") from e

        self.T = T
        self.P = P
        self.h = h
        self.d = d

        if phase == "twophase":
            Q = PropsSI_auto('Q', 'D', d, 'H', h, self.fluid)  # 0-1
            # saturated liquid and vapor properties at P
            h_l = PropsSI_auto('H', 'P', P, 'Q', 0, self.fluid)
            h_v = PropsSI_auto('H', 'P', P, 'Q', 1, self.fluid)
            d_l = PropsSI_auto('D', 'P', P, 'Q', 0, self.fluid)
            d_v = PropsSI_auto('D', 'P', P, 'Q', 1, self.fluid)

            self.Q = Q
            self.h_l = h_l
            self.h_v = h_v
            self.d_l = d_l
            self.d_v = d_v

            # masses
            self.m_v = Q * self.m
            self.m_l = self.m - self.m_v

            # fill level (volume fraction of liquid in tank)
            # liquid volume = m_l / d_l
            self.fill_level = (self.m_l / self.d_l) / self.V
        else:
            # single phase (liquid or gas)
            self.Q = None
            self.m_v = self.m if phase in ("gas", "supercritical") else 0.0
            self.m_l = self.m - self.m_v

            # set phase-specific properties equal to bulk
            self.h_l = self.h_v = self.h
            self.d_l = self.d_v = self.d
            self.fill_level = 1.0 if phase == "liquid" else 0.0

        # safe Cp/Cv/gamma/R in single-phase gas
        try:
            self.Cp = PropsSI_auto('CPMASS', 'D', self.d, 'H', self.h, self.fluid)
            self.Cv = PropsSI_auto('CVMASS', 'D', self.d, 'H', self.h, self.fluid)
            self.gamma = self.Cp / self.Cv if (self.Cv and self.Cp) else None
            self.R = self.Cp - self.Cv if (self.Cp and self.Cv) else None
        except Exception:
            self.Cp = self.Cv = self.gamma = self.R = None

    def update(self, mdot, Hdot, dt):
        """
        Updates node state based on an input mdot (kg/s), an input Hdot (J/s),
        as well as the sim timestep dt (s).
        """
        # apply updates
        self.m += mdot * dt
        self.H += Hdot * dt

        # numerical safety
        if self.m < 1e-12:
            self.m = 1e-12

        # recompute density and flash to get phase split
        d_new = self.m / self.V
        self._flash_from_DH(d_new, self.H)

        # debug print
        # print(f"{self.name}: t-update P={self.P:.1f} Pa, T={self.T:.3f} K, m={self.m:.6f} kg, m_l={self.m_l:.6f}, m_v={self.m_v:.6f}, Q={self.Q}")

    def log_state(self, t=0.0):
        """
        Log node state at each timestep throughout a network sim.
        """
        self.history["time"].append(t)
        self.history["Q"].append(self.Q)
        self.history["P"].append(self.P)
        self.history["T"].append(self.T)
        self.history["H"].append(self.H)
        self.history["h"].append(self.h)
        self.history["d"].append(self.d)
        self.history["m"].append(self.m)
        self.history["m_l"].append(self.m_l)
        self.history["m_v"].append(self.m_v)
        self.history["fill_level"].append(self.fill_level)


class Ambient(Node):
    """
    Subclass of Node to represnt ambient properties. 
    Unchanging regardless of updates into or out of it.
    """
    def __init__(self, fluid="Air", P=101325, T=293.15, name="ambient"):
        super().__init__(fluid, m=1.0, V=1.0, T=T, name=name)

        # Set fixed ambient conditions
        self.P = P
        self.T = T
        self.h = PropsSI_auto("H", "P", self.P, "T", self.T, fluid)
        self.H = self.h * self.m
        self.d = PropsSI_auto("D", "P", self.P, "T", self.T, fluid)

    def update(self, mdot, Hdot, dt):
        """
        Ignore mass/energy inflows, hold fixed at initial state.
        """        
        pass


class Manifold(Node):
    """
    Subclass of Node to represent a volumeless manifold.
    """
    # TODO
    pass


class Accumulator(Node):
    """
    Subclass of Node to represent an accumulator.
    """
    # TODO
    pass


class HeatExchanger(Node):
    """
    Subclass of Node to represnt a heat exchanger.
    """
    # TODO
    pass


class Chamber(Node):
    """
    Subclass of Node to represent a engine combustion chamber.
    """
    # TODO
    pass


class Connection():
    """
    Connection class. Defined by CdA (m^2), qdot (J/s), location on node (0-1), and state (open, closed).
    Initialized by CdA, qdot, location, and normal state.
    """
    def __init__(self, CdA, qdot=0.0, location=0.0, normal_state=True):
        self.CdA = CdA
        self.qdot = qdot
        self.location = location  # normalized height 0-1
        self.state = normal_state


    def mdot_Hdot(self, node1, node2):
        """
        Return mdot (kg/s), Hdot (J/s) where positive means mass/enthalpy flows node1 -> node2.
        Stream enthalpy is taken from donor's phase at the connection 'location'.
        """
        # check if connection is open
        if not self.state:
            return 0.0, 0.0
    
        dP = node1.P - node2.P
        if abs(dP) < 1e-12:
            return 0.0, 0.0

        # determine donor and receiver
        if dP > 0:
            donor, receiver = node1, node2
        else:
            donor, receiver = node2, node1

        # choose stream phase based on donor.fill_level and connection location
        # if fill_level > location -> stream draws from liquid; otherwise vapor
        if donor.fill_level > self.location:
            # stream is liquid from donor
            h_stream = donor.h_l
            # at receiver pressure
            d_stream = donor.d_l
        else:
            # stream is vapor from donor
            h_stream = donor.h_v
            d_stream = donor.d_v

        abs_dP = abs(dP)

        # allow choked formula only if donor is gas-like (no two-phase)
        donor_phase = CP.PhaseSI('D', donor.d, 'H', donor.h, donor.fluid)
        if donor_phase in ("gas", "supercritical") and donor.Cp and donor.Cv and donor.R:
            # compressible/choked flow estimate
            gamma = donor.gamma
            R = donor.R
            Tdon = donor.T
            crit_factor = ((gamma + 1.0) / 2.0) ** ( - (gamma + 1.0) / (2.0 * (gamma - 1.0)) )
            mdot_mag = self.CdA * donor.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor
        else:
            # incompressible-like orifice
            mdot_mag = self.CdA * np.sqrt(2.0 * max(d_stream, 1e-6) * abs_dP)
            Q_fill = PropsSI_auto('Q', 'P', receiver.P, 'H', h_stream, donor.fluid)
            # If inside 0–1, it's two-phase
            if 0 <= Q_fill <= 1:
                h_liq = PropsSI_auto('H', 'P', receiver.P, 'Q', 0, donor.fluid)
                h_vap = PropsSI_auto('H', 'P', receiver.P, 'Q', 1, donor.fluid)
                m_vap = mdot_mag * Q_fill
                m_liq = mdot_mag * (1 - Q_fill)
                # Hdot adjusted for flashing
                Hdot = m_vap * h_vap + m_liq * h_liq
            else:
                # single-phase at receiver P
                Hdot = mdot_mag * h_stream

        # sign convention: positive mdot -> node1 -> node2
        if donor is node1:
            mdot = mdot_mag
        else:
            mdot = - mdot_mag

        Hdot = mdot * h_stream + self.qdot
        return mdot, Hdot

    def mdot_Hdot_Dyer(self, node1, node2, choke_samples=24):
        """
        Compute mdot (kg/s) and Hdot (J/s) using a Dyer-style flashing relationship with automatic choking detection.
        Positive means mass/enthalpy flows node1 -> node2.

        choke_samples: number of downstream pressure samples to probe when searching for choked limit.
        """
        # closed?
        if not self.state:
            return 0.0, 0.0

        dP = node1.P - node2.P
        if abs(dP) < 1e-12:
            return 0.0, 0.0

        # donor/receiver selection (donor = higher pressure side)
        donor, receiver = (node1, node2) if dP > 0 else (node2, node1)
        Pu = donor.P
        Pd_actual = receiver.P
        abs_dP = abs(dP)

        # choose donor stream phase (based on fill height)
        if donor.fill_level > self.location:
            # drawing liquid from donor
            h_stream = donor.h_l
            rho_l = donor.d_l
            rho_v = donor.d_v
        else:
            # drawing vapor from donor
            h_stream = donor.h_v
            rho_l = donor.d_l
            rho_v = donor.d_v

        # donor bulk-phase classification for branching
        donor_phase = CP.PhaseSI('D', donor.d, 'H', donor.h, donor.fluid)

        # helper: compute Dyer-like mdot magnitude for a given downstream pressure Pd
        def dyer_mdot_for_Pd(Pd):
            # if Pd >= Pu, no forward flow
            if Pd >= Pu:
                return 0.0

            # fetch sat props at upstream pressure (Pu)
            try:
                h_lsat_up = PropsSI_auto('H', 'P', Pu, 'Q', 0, donor.fluid)
                h_vsat_up = PropsSI_auto('H', 'P', Pu, 'Q', 1, donor.fluid)
                h_fg_up = max(h_vsat_up - h_lsat_up, 1e-12)
            except Exception:
                # if CoolProp cannot provide saturated props at Pu, fallback to simple orifice
                return self.CdA * np.sqrt(2.0 * max(rho_l, 1e-6) * (Pu - Pd))

            # donor specific enthalpy (bulk)
            hu = donor.h

            # available flashing fraction based on upstream energy
            x_available = (hu - h_lsat_up) / h_fg_up
            x_available = float(np.clip(x_available, 0.0, 1.0))

            # Dyer-style empirical factor to account for two-phase expansion:
            # build a blending based on density ratio and available quality.
            Rho_ratio_sqrt = np.sqrt(max(rho_l / max(rho_v, 1e-12), 1e-12))
            Y = (1.0 + Rho_ratio_sqrt * x_available) / (1.0 + Rho_ratio_sqrt)

            # base mass-flux-like term
            deltaP = max(Pu - Pd, 0.0)
            G = np.sqrt(2.0 * rho_l * deltaP * Y)  # [kg/(m^2 s)] approximate

            mdot_mag = self.CdA * G
            return float(mdot_mag)

        # --- compute candidate mdot for actual Pd (no choke cap yet) ---
        if donor_phase in ("gas", "supercritical"):
            # single-phase gas: use choked-gas expression (isentropic approximate)
            if donor.Cp and donor.Cv and donor.R and donor.gamma:
                gamma = donor.gamma
                R = donor.R
                Tdon = donor.T
                crit_factor = ((gamma + 1.0) / 2.0) ** ( - (gamma + 1.0) / (2.0 * (gamma - 1.0)) )
                mdot_mag_nominal = self.CdA * donor.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor
            else:
                mdot_mag_nominal = self.CdA * np.sqrt(2.0 * max((rho_l if rho_l is not None else 1.0), 1e-6) * abs_dP)
        else:
            # flashing/two-phase or liquid: use Dyer-like mdot at actual Pd
            mdot_mag_nominal = dyer_mdot_for_Pd(Pd_actual)

        # --- choking detection & cap (only meaningful for two-phase/flashing donors) ---
        mdot_mag = mdot_mag_nominal
        if donor_phase in ("twophase", "liquid", "supercritical_liquid"):
            # scan Pd values from Pd_actual down to a small fraction of Pu to find max possible mdot
            Pd_min = max(Pu * 1e-4, 1.0)  # avoid zero
            # log-spaced samples between Pd_actual and Pd_min (including Pd_actual)
            # ensure Pd_actual >= Pd_min
            low = min(Pd_actual, Pd_min)
            high = max(Pd_actual, Pd_min)
            # sample logspace between high and low (descending)
            Pd_samples = np.unique(
                np.concatenate((
                    np.linspace(Pd_actual, Pu * 1e-4, choke_samples),
                    np.geomspace(max(Pd_actual, 1e-3), max(Pu * 1e-4, 1e-3), choke_samples)
                ))
            )
            # compute mdot for each Pd sample
            mdot_vals = np.array([dyer_mdot_for_Pd(Pd) for Pd in Pd_samples])
            mdot_max = float(np.max(mdot_vals)) if mdot_vals.size else mdot_mag_nominal
            # if the nominal mdot exceeds the maximum achievable massflux, cap it
            if mdot_mag_nominal > mdot_max:
                mdot_mag = mdot_max
            else:
                mdot_mag = mdot_mag_nominal

        # --- compute Hdot consistent with receiver-phase flashing at Pd_actual ---
        # First, find stream enthalpy at donor (we already have h_stream as donor.h_l or h_v)
        # Now evaluate quality at receiver pressure for that h_stream
        Hdot = 0.0
        try:
            Q_recv = PropsSI_auto('Q', 'P', Pd_actual, 'H', h_stream, donor.fluid)
        except Exception:
            # If CoolProp can't compute Q (e.g. single phase), fallback to single-phase convective enthalpy
            Q_recv = None

        if Q_recv is not None and (0.0 <= Q_recv <= 1.0):
            # Stream flashed at receiver pressure: split enthalpy into saturated contributions
            h_liq_recv = PropsSI_auto('H', 'P', Pd_actual, 'Q', 0, donor.fluid)
            h_vap_recv = PropsSI_auto('H', 'P', Pd_actual, 'Q', 1, donor.fluid)
            m_vap = mdot_mag * Q_recv
            m_liq = mdot_mag * (1.0 - Q_recv)
            Hdot = m_vap * h_vap_recv + m_liq * h_liq_recv
        else:
            # single-phase stream at receiver P: use donor stream enthalpy as convected
            Hdot = mdot_mag * h_stream

        # add heat transfer through the connection
        Hdot = Hdot + self.qdot

        # sign convention: positive mdot -> node1 -> node2
        if donor is node1:
            mdot = mdot_mag
        else:
            mdot = - mdot_mag

        return mdot, Hdot


class Regulator(Connection):
    def __init__(self, CdA, set_pressure, droop_curve=None, qdot=0.0, location=0.0, normal_state=True):
        """
        A pressure regulator connection that limits downstream pressure. Defined by: CdA (m^2), set_pressure (Pa),
        droop_curve (function that maps mdot -> pressure drop (Pa)), qdot (J/s), location (0-1), state (open, closed).
        """
        super().__init__(CdA, qdot, location, normal_state)
        self.set_pressure = set_pressure
        self.droop_curve = droop_curve  # function handle: f(mdot) -> ΔP droop

    def mdot_Hdot(self, node1, node2):
        """
        Computes mdot and Hdot across the regulator.
        The regulator limits downstream pressure to set_pressure (minus droop if defined).
        Args:
            node1, node1 (Node): nodes connected by this connection
        """
        if not self.state:
            return 0.0, 0.0

        # Determine upstream and downstream
        dP = node1.P - node2.P
        if abs(dP) < 1e-12:
            return 0.0, 0.0

        if dP > 0:
            upstream, downstream = node1, node2
        else:
            upstream, downstream = node2, node1

        # Target downstream pressure
        P_down_target = self.set_pressure

        # Apply droop curve if defined
        if self.droop_curve is not None:
            # iterative droop correction: assume mdot ≈ previous mdot, or start with 0
            # droop curve returns positive ΔP loss at higher flows
            P_down_target -= self.droop_curve(abs(dP))  

        # Clamp downstream pressure to not exceed target
        if downstream.P < P_down_target:
            # regulator closed: no flow (receiver pressure too low)
            return 0.0, 0.0
        else:
            # regulator open: limit flow so that downstream ≈ setpoint
            effective_dP = max(upstream.P - P_down_target, 0.0)

        # Now use inherited orifice logic for the flow
        if upstream.fill_level > self.location:
            h_stream = upstream.h_l
            d_stream = upstream.d_l
        else:
            h_stream = upstream.h_v
            d_stream = upstream.d_v

        donor_phase = CP.PhaseSI('D', upstream.d, 'H', upstream.h, upstream.fluid)
        if donor_phase in ("gas", "supercritical") and upstream.Cp and upstream.Cv and upstream.R:
            gamma = upstream.gamma
            R = upstream.R
            Tdon = upstream.T
            crit_factor = ((gamma + 1.0) / 2.0) ** ( - (gamma + 1.0) / (2.0 * (gamma - 1.0)) )
            mdot_mag = self.CdA * upstream.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor
        else:
            mdot_mag = self.CdA * np.sqrt(2.0 * max(d_stream, 1e-6) * effective_dP)

        # Sign convention: positive mdot -> node1 -> node2
        mdot = mdot_mag if upstream is node1 else -mdot_mag

        Hdot = mdot * h_stream + self.qdot
        return mdot, Hdot


class Valve(Connection):
    """
    Subclass of Connection to represent a valve.
    """
    # TODO
    pass


class SharpEdgedOrifice(Connection):
    """
    Subclass of Node to represent a sharp-edged orifice.
    """
    # TODO
    pass


class Engine(Connection):
    """
    Subclass of Connection to represent an engine.
    """    
    # TODO
    pass


class Injector(Connection):
    """
    Subclass of Connection to represent an injector.
    """
    # TODO
    pass


class ThrottleValve(Connection):
    """
    Subclass of Connection to represent a throttle valve.
    """
    # TODO
    pass


class CheckValve(Connection):
    """
    Subclass of Connection to represent a check valve.
    """
    # TODO
    pass


class Pump(Connection):
    """
    Subclass of Connection to represent a pump.
    """
    # TODO
    pass


class Network():
    """
    Network class. Defined by a graph of connections and nodes representing a fluid network.
    Syntax is as follows: {connection1: (node1, node2), ...}
    """
    def __init__(self, graph):
        self.graph = graph  # {connection: (node1, node2)}

    def sim(self, t, dt, actions={}, verbose_steps=5):
        """
        Runs a transient sim over a fluid network. Takes in total sim runtime, timestep, and an action profile of
        valve set-states and timings.
        Args:
            t: total sim runtime [s] (float)
            dt: timestep duration [s] (float)
            actions: {time [s]: (connection, set_state), ...}, action time must be of the same order of precision as dt.
        """
        steps = int(t / dt)
        for i in range(steps):
            time_now = i * dt
            if time_now in actions:
                actions[time_now][0].state = actions[time_now][1]

            # compute all flows first
            mdot_contrib = {node: 0.0 for _, pair in self.graph.items() for node in pair}
            Hdot_contrib = {node: 0.0 for _, pair in self.graph.items() for node in pair}

            for conn, (n1, n2) in self.graph.items():
                mdot, Hdot = conn.mdot_Hdot(n1, n2)
                # mdot positive => n1 -> n2
                mdot_contrib[n1] -= mdot
                mdot_contrib[n2] += mdot
                Hdot_contrib[n1] -= Hdot
                Hdot_contrib[n2] += Hdot

            # update nodes simultaneously
            for node in list(mdot_contrib.keys()):
                node.update(mdot_contrib[node], Hdot_contrib[node], dt)
                if i < verbose_steps:
                # print only first few steps to avoid log overload
                    print(f"[t={time_now:.4f}] {node.name} mdot_net={mdot_contrib[node]:.6f}, Hdot_net={Hdot_contrib[node]:.3f}")
                node.log_state(time_now)

    def plot_nodes(self, nodes, units="SI"):
        """
        Plots pressure, temperature, mass, density, quality, and fill level vs time
        for each node as subplots in one figure per node.

        Args:
            time [s] (array): array of time values
            nodes (list): list of Node objects with histories recorded
            units (str): SI or E
        """
        if units == "E":
            psi = 1 / 6894.75729
        else:
            psi = 1
        for node in nodes:
            time = node.history['time']
            fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
            axs = axs.flatten()  # make it a flat list for easy indexing
            fig.suptitle(f"Node: {node.name}", fontsize=14)

            # Pressure
            axs[0].plot(time, node.history['P']*psi)
            axs[0].set_ylabel("Pressure [Pa]")

            # Temperature
            axs[1].plot(time, node.history['T'])
            axs[1].set_ylabel("Temperature [K]")

            # Mass
            axs[2].plot(time, node.history['m'])
            axs[2].set_ylabel("Mass [kg]")

            # Density
            axs[3].plot(time, node.history['d'])
            axs[3].set_ylabel("Density [kg/m³]")

            # Quality
            axs[4].plot(time, node.history['Q'])
            axs[4].set_ylabel("Quality [-]")
            axs[4].set_ylim(0, 1)  # quality is between 0–1

            # Fill level
            axs[5].plot(time, node.history['fill_level'])
            axs[5].set_ylabel("Fill level [m]")
            axs[5].set_xlabel("Time [s]")

            plt.tight_layout(rect=[0, 0, 1, 0.96])
            plt.show()

    def plot_nodes_overlay(self, nodes, title="Node Comparison", units="SI"):
        """
        Overlay plots of pressure, temperature, mass, density, quality,
        and fill level vs time for all nodes on the same set of subplots.
        Args:
            nodes (list): list of Node objects with histories recorded
            title (str): plot title
            units (str): SI or E
        """
        fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
        axs = axs.flatten()
        fig.suptitle(title, fontsize=14)

        # Loop over nodes and add each to the plots
        for node in nodes:
            time = node.history['time']
            if units == "E":
                axs[0].plot(time, np.array(node.history['P']) / 6894.75729, label=node.name)
                axs[1].plot(time, (np.array(node.history['T']) - 273.15) * 1.8 + 32, label=node.name)
            else:
                axs[0].plot(time, node.history['P'], label=node.name)
                axs[1].plot(time, node.history['T'], label=node.name)
            axs[2].plot(time, node.history['m'], label=node.name)
            axs[3].plot(time, node.history['d'], label=node.name)
            axs[4].plot(time, node.history['Q'], label=node.name)
            axs[5].plot(time, node.history['fill_level'], label=node.name)

        # Labels
        if units == "E":
            axs[0].set_ylabel("Pressure [psi]")
            axs[1].set_ylabel("Temperature [F]")
        else:
            axs[0].set_ylabel("Pressure [Pa]")
            axs[1].set_ylabel("Temperature [K]")
        axs[2].set_ylabel("Mass [kg]")
        axs[3].set_ylabel("Density [kg/m³]")
        axs[4].set_ylabel("Quality [-]")
        axs[4].set_ylim(0, 1)
        axs[5].set_ylabel("Fill level [-]")
        axs[5].set_xlabel("Time [s]")

        # Add legends
        for ax in axs:
            ax.legend()
            ax.grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()