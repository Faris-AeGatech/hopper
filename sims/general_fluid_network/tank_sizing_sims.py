from general_fluid_network import Node, Ambient, Connection, Regulator, Network
#################### TEST CONFIGS #########################

# GAS PULL
# test_node1 = Node("N2O", 35, 60, 293, "Liquid Pull")
# test_node2 = Node("N2O", 35, 60, 293, "Gas Pull")
# amb_node = Ambient()
# test_connection1 = Connection(0.0000063, 0, 0)
# test_connection2 = Connection(0.00001521, 0, 1)
# test_connection3 = Connection(0.0000001742, 0, 0)
# test_network = Network({test_connection1: (test_node1, amb_node), test_connection2: (test_node2, amb_node)})
# test_network.sim(60, 1)
# test_network.plot_nodes_overlay((test_node1, test_node2), units="E")

# LIQUID PULL FILL
# target mass: 36 kg
n = 5 # num bottles
# vehicle_tank = Node("N2O", 1, 40, 293, "Vehcle Tank")
# fill_tanks = Node("N2O", 9.0718 * n, 13.4 * n, 293, "Fill Tanks")
# amb_node = Ambient()
# fill_line = Connection(0.000006, 0, 0)
# vent_line = Connection(0.0000003, 0, 1, False)
# omv = Connection(0.000008, 0, 0, False)
# test_network = Network({fill_line: (fill_tanks, vehicle_tank), vent_line: (vehicle_tank, amb_node), omv: (vehicle_tank, amb_node)})
# test_network.sim(940, 1, {250: (vent_line, True), 870: (fill_line, False), 880: (omv, True)})
# test_network.plot_nodes_overlay((fill_tanks, vehicle_tank), title="LIQUID PULL 40L", units="E")101325

# LIQUID PULL
vehicle_tank = Node("N2O", 36, 60, 293, "Liquid Pull")
amb_node = Ambient(P=101325*350/14.7)
fluid_system = Connection(0.000009, 0, 0)
test_network = Network({fluid_system: (vehicle_tank, amb_node)})
test_network.sim(60, 1)
test_network.plot_nodes_overlay((vehicle_tank, vehicle_tank), units="E")

# LIQUID PULL SMALL TANK
# n = 1
# vehicle_tank = Node("N2O", 36, 50, 293, "Vehicle Tank")
# vt2 = Node("N2O", 9.0718 * n, 13.4 * n, 293, "Bottle")
# amb_node = Ambient(P=101325*300/14.7)
# fluid_system = Connection(0.000009, 0, 0)
# fluid_system2 = Connection(0.000009, 0, 0)
# test_network = Network({fluid_system: (vehicle_tank, amb_node), fluid_system2: (vt2, amb_node)})
# test_network.sim(80, 1)
# test_network.plot_nodes_overlay((vehicle_tank, vt2), units="E")


# SUBCOOLED COPV 
# copv = Node("Nitrogen", 2.4, 6.61, 293, "COPV")
# vehicle_tank = Node("N2O", 36, 60, 287, "Liquid Pull")
# amb_node = Ambient()
# tank_substitute = Ambient(fluid="Nitrogen", P=101325*400/14.7)
# reg = Regulator(0.00000063, 101325*400/14.7)
# fluid_system = Connection(0.0000063, 0, 0)
# network = Network({reg: (copv, vehicle_tank), fluid_system: (vehicle_tank, amb_node)})
# network.sim(60, 1)
# network.plot_nodes_overlay((copv, vehicle_tank), units="E")

# DARCY SPACE FILL VALIDATION
# target mass: 108.86 kg
# n = 24 # num bottles
# vehicle_tank = Node("N2O", 1, 106.5, 293, "Vehicle Tank")
# fill_tanks = Node("N2O", 9.0718 * n, 13.4 * n, 293, "Fill Tanks")
# amb_node = Ambient()
# fill_line = Connection(0.000006, 0, 0)
# vent_line = Connection(0.0000003, 0, 1, False)
# omv = Connection(0.000008, 0, 0, False)
# test_network = Network({fill_line: (fill_tanks, vehicle_tank), vent_line: (vehicle_tank, amb_node), omv: (vehicle_tank, amb_node)})
# test_network.sim(1000, 1, {400: (vent_line, True)})
# test_network.plot_nodes_overlay((fill_tanks, vehicle_tank), title="Darcy Space Validation", units="E")

