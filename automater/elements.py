class GenericClass(object):
	pass

class Element:
	"""Atomic information from periodic table"""
	def __init__(self, symbol = False):
		if not symbol:
			exit('New elements in elements.py must be instantiated with an elemental symbol')
		self.symbol = symbol
		self.longName = False
		self.mass = False
		self.pseudoList = False
		self.longName = False
		self.number = False
		self.mass = False
		self.radius = False # atomic radius in picometer
		self.electonegativity = False
		self.series = False
		self.group = GenericClass()
		self.group.number = False
		self.group.CAS = False
		self.group.name = False
		self.period = False
		self.valence = GenericClass()
		self.valence.n = False
		self.valence.shellNumber = False
		self.valence.l = False
		self.valence.shell = False
		self.valence.subshell = False
		self.valence.orbital = False
		self.valence.configuration = False
		self.valence.electrons = False
		self.structures = GenericClass()
		self.structures.primitive = False
		self.structures.conventional = False

H = h = Hydrogen = hydrogen = Element('H') # Instantiate element with atomic symbol 
H.longName = 'Hydrogen' # Full English name of the element
H.number = 1 # Atomic number 
H.mass = 1.008 # Mass in amu
H.radius = 53 # atomic radius in picometer
H.electonegativity = 2.20 
H.series = 'Nonmetal' 
H.group = GenericClass()
H.group.number = 1
H.group.CAS = 'IA'
H.group.name = 'Alkali Metal'
H.period = 1
H.valence =  GenericClass()
H.valence.n = H.valence.shellNumber = 1
H.valence.l = 0
H.valence.shell = 'K'
H.valence.subshell = 's'
H.valence.orbital = '1s'
H.valence.configuration = '1s1'
H.valence.electrons = 1
H.structures = GenericClass()
H.structures.primitive = 'H2_mp-632250_primitive.cif'
H.structures.conventional = 'H2_mp-632250_conventional_standard.cif'

He = he = Helium = helium = Element('He')
He.longName = 'Helium' 
He.number = 2
He.mass = 4.002602
He.radius =  31 # atomic radius in picometer
He.electonegativity = False
He.series = 'Noble Gas'
He.group = GenericClass()
He.group.number = 18
He.group.CAS = 'VIIIA'
He.group.name = 'Noble Gas'
He.period = 1
He.valence =  GenericClass()
He.valence.n = H.valence.shellNumber = 1
He.valence.l = 0
He.valence.shell = 'K'
He.valence.subshell = 's'
He.valence.orbital = '1s'
He.valence.configuration = '1s2'
He.valence.electrons = 2
He.structures = GenericClass()
He.structures.primitive = 'He_mp-23158_primitive.cif'
He.structures.conventional = 'He_mp-23158_conventional_standard.cif'

Li = li = Lithium = lithium = Element('Li')
Li.longName = 'Lithium'
Li.number = 3
Li.mass = 6.94
Li.radius = 167 # atomic radius in picometer
Li.electonegativity = 0.98
Li.series = 'Alkali Metal'
Li.group = GenericClass()
Li.group.number = 1
Li.group.CAS = 'IA'
Li.group.name = 'Alkali Metal'
Li.period = 2
Li.valence =  GenericClass()
Li.valence.n = H.valence.shellNumber = 2
Li.valence.l = 0
Li.valence.shell = 'L'
Li.valence.subshell = 's'
Li.valence.orbital = '2s'
Li.valence.configuration = '2s1'
Li.valence.electrons = 1
Li.structures = GenericClass()
Li.structures.primitive = 'Li_mp-135_primitive.cif'
Li.structures.conventional = 'Li_mp-135_conventional_standard.cif'


As = Arsenic = arsenic = Element('As') # Not using "as" as an alias, because it's a keyword in python
As.longName = 'Arsenic'
As.number = 33
As.mass = 74.92160
As.radius =  114 # atomic radius in picometer
As.electonegativity = 2.18
As.series = 'Metalloid'
As.group = GenericClass()
As.group.number = 15
As.group.CAS = 'VA'
As.group.name = 'Pnictogen'
As.period = 4
As.valence =  GenericClass()
As.valence.n = H.valence.shellNumber = 4
As.valence.l = 1
As.valence.shell = 'N'
As.valence.subshell = 'p'
As.valence.orbital = '4p'
As.valence.configuration = '3d10-4s2-4p3'
As.valence.electrons = 5
As.structures = GenericClass()
As.structures.primitive = 'As_mp-11_primitive.cif'
As.structures.conventional = 'As_mp-11_conventional_standard.cif'

Cd = cd = Cadmium = cadmium = Element('Cd')
Cd.longName = 'Cadmium'
Cd.number = 48
Cd.mass = 112.411
Cd.radius =  161 # atomic radius in picometer
Cd.electonegativity = 1.69
Cd.series = 'Transition Metal'
Cd.group = GenericClass()
Cd.group.number = 12
Cd.group.CAS = 'IIB'
Cd.group.name = 'Volatile Metal'
Cd.period = 5
Cd.valence =  GenericClass()
Cd.valence.n = H.valence.shellNumber = 5
Cd.valence.l = 0
Cd.valence.shell = 'O'
Cd.valence.subshell = 's'
Cd.valence.orbital = '5s'
Cd.valence.configuration = '4d10-5s2'
Cd.valence.electrons = 12
Cd.structures = GenericClass()
Cd.structures.primitive = 'Cd_mp-94_primitive.cif'
Cd.structures.conventional = 'Cd_mp-94_conventional_standard.cif'

Te = te = Tellurium = tellurium = Element('Te')
Te.longName = 'Tellurium'
Te.number = 52
Te.mass = 127.60
Te.radius = 123 # atomic radius in picometer
Te.electonegativity = 2.1
Te.series = 'Metalloid'
Te.group = GenericClass()
Te.group.number = 16
Te.group.CAS = 'VIA'
Te.group.name = 'Chalcogen'
Te.period = 5
Te.valence =  GenericClass()
Te.valence.n = H.valence.shellNumber = 5
Te.valence.l = 1
Te.valence.shell = 'O'
Te.valence.subshell = 'p'
Te.valence.orbital = '5p'
Te.valence.configuration = '4d10-5s2-5p4'
Te.valence.electrons = 6
Te.structures = GenericClass()
Te.structures.primitive = 'Te_mp-19_primitive.cif'
Te.structures.conventional = 'Te_mp-19_conventional_standard.cif'

'''
= #Name = Element()
.longName = #
.number = 1
.mass = 1.008
.radius = 53 # atomic radius in picometer
.electonegativity = 2.20
.series = 'Nonmetal'
.group = GenericClass()
.group.number = 1
.group.CAS = 'IA'
.group.name = 'Alkali Metal'
.period = 1
.valence =  GenericClass()
.valence.n = H.valence.shellNumber = 1
.valence.l = 0
.valence.shell = 'K'
.valence.subshell = 's'
.valence.orbital = '1s'
.valence.configuration = '1s1'
.valence.electrons = '1'
.structures = GenericClass()
.structures.primitive = 'H2_mp-632250_primitive.cif'
.structures.conventional = 'H2_mp-632250_conventional_standard.cif'
'''



'''
1	H	Hydrogen
2	He	Helium
3	Li	Lithium
4	Be	Beryllium
5	B	Boron
6	C	Carbon
7	N	Nitrogen
8	O	Oxygen
9	F	Fluorine
10	Ne	Neon
11	Na	Sodium
12	Mg	Magnesium
13	Al	Aluminum
14	Si	Silicon
15	P	Phosphorus
16	S	Sulfur
17	Cl	Chlorine
18	Ar	Argon
19	K	Potassium
20	Ca	Calcium
21	Sc	Scandium
22	Ti	Titanium
23	V	Vanadium
24	Cr	Chromium
25	Mn	Manganese
26	Fe	Iron
27	Co	Cobalt
28	Ni	Nickel
29	Cu	Copper
30	Zn	Zinc
31	Ga	Gallium
32	Ge	Germanium
33	As	Arsenic
34	Se	Selenium
35	Br	Bromine
36	Kr	Krypton
37	Rb	Rubidium
38	Sr	Strontium
39	Y	Yttrium
40	Zr	Zirconium
41	Nb	Niobium
42	Mo	Molybdenum
43	Tc	Technetium
44	Ru	Ruthenium
45	Rh	Rhodium
46	Pd	Palladium
47	Ag	Silver
48	Cd	Cadmium
49	In	Indium
50	Sn	Tin
51	Sb	Antimony
52	Te	Tellurium
53	I	Iodine
54	Xe	Xenon
55	Cs	Cesium
56	Ba	Barium
57	La	Lanthanum
58	Ce	Cerium
59	Pr	Praseodymium
60	Nd	Neodymium
61	Pm	Promethium
62	Sm	Samarium
63	Eu	Europium
64	Gd	Gadolinium
65	Tb	Terbium
66	Dy	Dysprosium
67	Ho	Holmium
68	Er	Erbium
69	Tm	Thulium
70	Yb	Ytterbium
71	Lu	Lutetium
72	Hf	Hafnium
73	Ta	Tantalum
74	W	Tungsten
75	Re	Rhenium
76	Os	Osmium
77	Ir	Iridium
78	Pt	Platinum
79	Au	Gold
80	Hg	Mercury
81	Tl	Thallium
82	Pb	Lead
83	Bi	Bismuth
84	Po	Polonium
85	At	Astatine
86	Rn	Radon
87	Fr	Francium
88	Ra	Radium
89	Ac	Actinium
90	Th	Thorium
91	Pa	Protactinium
92	U	Uranium
93	Np	Neptunium
94	Pu	Plutonium
95	Am	Americium
96	Cm	Curium
97	Bk	Berkelium
98	Cf	Californium
99	Es	Einsteinium
100	Fm	Fermium
101	Md	Mendelevium
102	No	Nobelium
103	Lr	Lawrencium
104	Rf	Rutherfordium
105	Db	Dubnium
106	Sg	Seaborgium
107	Bh	Bohrium
108	Hs	Hassium
109	Mt	Meitnerium
110	Ds	Darmstadtium
111	Rg	Roentgenium
112	Cn	Copernicium
113	Nh	Nihonium
114	Fl	Flerovium
115	Mc	Moscovium
116	Lv	Livermorium
117	Ts	Tennessine
118	Og	Oganesson
'''
'''
PSEUDOPOTENTIALS
Si.blyp-hgh.UPF
Si.pz-vbc.UPF
Si.pw91-n-van.UPF
Cd.pw-mt_fhi.UPF
Cd.pbe-mt_fhi.UPF
Cd.pbe-dn-rrkjus_psl.0.2.UPF
Cd.pbe-dn-kjpaw_psl.0.2.UPF
Cd.pbesol-dn-rrkjus_psl.0.2.UPF
Cd.pbesol-dn-kjpaw_psl.0.2.UPF
Te.pw-mt_fhi.UPF
Te.pbe-mt_fhi.UPF
Te.pbesol-dn-kjpaw_psl.0.2.2.UPF
Te.pbe-dn-rrkjus_psl.0.2.2.UPF
Te.pbe-dn-kjpaw_psl.0.2.2.UPF
Te.pbesol-dn-rrkjus_psl.0.2.2.UPF
As.pbesol-n-kjpaw_psl.0.2.UPF # PAW PBESOL
As.pbe-n-kjpaw_psl.0.2.UPF
As.pbe-mt_fhi.UPF 
Cd.pbe-dn-kjpaw_psl.0.2.UPF
'''