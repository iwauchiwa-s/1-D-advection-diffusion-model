[header] for 32O2, temp read, TDF none
2 ! isys = sim. type (1:boundary atm., 2:stratosphere)
100 ! nzd = number of layer
1000.0 ! dzd1 = layer thickness [m]
10.0 ! dt = time step [s]
1901 2000 12 31 ! iys, iyp, imn, idays = start year, stop year/month/day
32.0e-3 1.268 0.00 ! amsg, sfd, tdf = molar mass[kg/mol], diff. scale factor, thermal diff. factor of gas
0 0.208471722 ! iatc, atcc = flag for upper boundary scenario & constant value
0   ! irinit = flag of initial read from file
1 0.4 ! irdedy, dedy0 = flag of Dedy read from file[m2/s]
1 273.0 ! irtemp, temp0 = flag of temp read from file[K]
0 0.0 1.0e-5 ! irwind, w0 = flag of w read from file[m/s]

for 32O2
32.0e-3 1.268 0.00 ! amsg, sfd, tdf = molar mass[kg/mol], diff. scale factor, thermal diff. factor of gas
0 0.208471722 ! iatc, atcc = flag for upper boundary scenario & constant value

for 34O2
34.0e-3 1.25160 0.00 ! amsg, sfd, tdf = molar mass[kg/mol], diff. scale factor, thermal diff. factor of gas 
0 0.000840072 ! iatc, atcc = flag for upper boundary scenario & constant value


note
natural existance ratio
16O=0.9976
17O=0.00039
18O=0.00201
->
32O2=0.99520576
34O2=18O-16O=0.004010352

atm. O2 total conc. = 0.209476

32O2 conc.=0.208471722
34O2 conc.=0.000840072
as zero permil
