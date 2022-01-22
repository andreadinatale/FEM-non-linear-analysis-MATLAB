function INPUT = inputNL_strip_hole_120_els
%
% Structure INPUT
%     integration_pts: number of integration points
%            elements: nodal connectivity
%               nodes: nodal coordinates
%                   E: Young's modulus
%                  nu: Poisson's ratio
%                   t: thickness
%            mat_type: material type (name from library)
%                      'neohookean': incompressible NH model; only for nonlinear analysis 
%                      'plane_stress': linear plane stress law model 
%                      'plane_stress_linear': linear plane stress law for *linear* analysis 
%                load: [node_id component magnitude]
%                 spc: [node_id component]
% ... additional info for nonlinear analysis
% 
%            sol_type: solution procedure
%                      'linear'
%                      'nonlinear'
%          lambda_max: max load parameter (only for nonlinear analysis; any number is ok for linear analysis)
%             dlambda: load increment  (only for nonlinear analysis; any number is ok for linear analysis)
%        norm_res_max: tolerance for NR iterations (only for nonlinear analysis; any number is ok for linear analysis)
%            nitermax: max number for NR iterations (only for nonlinear analysis; any number is ok for linear analysis)

% -- Init
INPUT = struct();

% -- Elements
INPUT.integration_pts = 2;
% INPUT.integration_barlow = 1; %optional field: to recover stresses with Barlow set to 1; otherwise omit or set to 0;

INPUT.elements = [   1     2     8     7
  2     3     9     8
  3     4    10     9
  4     5    11    10
  5     6    12    11
  7     8    14    13
  8     9    15    14
  9    10    16    15
 10    11    17    16
 11    12    18    17
 13    14    20    19
 14    15    21    20
 15    16    22    21
 16    17    23    22
 17    18    24    23
 19    20    26    25
 20    21    27    26
 21    22    28    27
 22    23    29    28
 23    24    30    29
 25    26    32    31
 26    27    33    32
 27    28    34    33
 28    29    35    34
 29    30    36    35
 31    32    38    37
 32    33    39    38
 33    34    40    39
 34    35    41    40
 35    36    42    41
 37    38    44    43
 38    39    45    44
 39    40    46    45
 40    41    47    46
 41    42    48    47
 43    44    50    49
 44    45    51    50
 45    46    52    51
 46    47    53    52
 47    48    54    53
 49    50    56    55
 50    51    57    56
 51    52    58    57
 52    53    59    58
 53    54    60    59
 55    56    62    61
 56    57    63    62
 57    58    64    63
 58    59    65    64
 59    60    66    65
 67    68    79    78
 68    69    80    79
 69    70    81    80
 70    71    82    81
 71    72    83    82
 72    73    84    83
 73    74    85    84
 74    75    86    85
 75    76    87    86
 76    77    88    87
 78    79    90    89
 79    80    91    90
 80    81    92    91
 81    82    93    92
 82    83    94    93
 83    84    95    94
 84    85    96    95
 85    86    97    96
 86    87    98    97
 87    88    99    98
 89    90   101   100
 90    91   102   101
 91    92   103   102
 92    93   104   103
 93    94   105   104
 94    95   106   105
 95    96   107   106
 96    97   108   107
 97    98   109   108
 98    99   110   109
100   101   112   111
101   102   113   112
102   103   114   113
103   104   115   114
104   105   116   115
105   106   117   116
106   107   118   117
107   108   119   118
108   109   120   119
109   110   121   120
111   112    12     6
112   113    18    12
113   114    24    18
114   115    30    24
115   116    36    30
116   117    42    36
117   118    48    42
118   119    54    48
119   120    60    54
120   121    66    60 ];

% -- Nodes
INPUT.nodes = [  0.0000     0.2500
0.0391     0.2469
0.0772     0.2378
0.1135     0.2227
0.1470     0.2022
0.1768     0.1768
0.0000     0.5500
0.1002     0.5472
0.1995     0.5390
0.2972     0.5255
0.3923     0.5070
0.4841     0.4841
0.0000     0.8500
0.1613     0.8475
0.3218     0.8402
0.4808     0.8282
0.6376     0.8118
0.7914     0.7914
0.0000     1.1500
0.2224     1.1478
0.4441     1.1414
0.6645     1.1309
0.8829     1.1166
1.0987     1.0987
0.0000     1.4500
0.2835     1.4482
0.5663     1.4427
0.8481     1.4336
1.1282     1.4213
1.4061     1.4061
0.0000     1.7500
0.3445     1.7485
0.6886     1.7439
1.0318     1.7364
1.3735     1.7261
1.7134     1.7134
0.0000     2.0500
0.4056     2.0488
0.8109     2.0451
1.2154     2.0391
1.6188     2.0309
2.0207     2.0207
0.0000     2.3500
0.4667     2.3491
0.9332     2.3463
1.3991     2.3418
1.8641     2.3357
2.3280     2.3280
0.0000     2.6500
0.5278     2.6494
1.0554     2.6476
1.5827     2.6446
2.1094     2.6404
2.6354     2.6354
0.0000     2.9500
0.5889     2.9497
1.1777     2.9488
1.7664     2.9473
2.3547     2.9452
2.9427     2.9427
0.0000     3.2500
0.6500     3.2500
1.3000     3.2500
1.9500     3.2500
2.6000     3.2500
3.2500     3.2500
0.2500     0.0000
0.5500     0.0000
0.8500     0.0000
1.1500     0.0000
1.4500     0.0000
1.7500     0.0000
2.0500     0.0000
2.3500     0.0000
2.6500     0.0000
2.9500     0.0000
3.2500     0.0000
0.2469     0.0391
0.5472     0.1002
0.8475     0.1613
1.1478     0.2224
1.4482     0.2835
1.7485     0.3445
2.0488     0.4056
2.3491     0.4667
2.6494     0.5278
2.9497     0.5889
3.2500     0.6500
0.2378     0.0772
0.5390     0.1995
0.8402     0.3218
1.1414     0.4441
1.4427     0.5663
1.7439     0.6886
2.0451     0.8109
2.3463     0.9332
2.6476     1.0554
2.9488     1.1777
3.2500     1.3000
0.2227     0.1135
0.5255     0.2972
0.8282     0.4808
1.1309     0.6645
1.4336     0.8481
1.7364     1.0318
2.0391     1.2154
2.3418     1.3991
2.6445     1.5827
2.9473     1.7664
3.2500     1.9500
0.2022     0.1470
0.5070     0.3923
0.8118     0.6376
1.1166     0.8829
1.4213     1.1282
1.7261     1.3735
2.0309     1.6188
2.3357     1.8641
2.6404     2.1094
2.9452     2.3547
3.2500     2.6000 ];

% --- Material properties
INPUT.mat_type = 'neohookean'; 
INPUT.E  = 1.2675;
INPUT.nu = .5;
INPUT.t  = 0.079;
  
% -- Loading conditions
Fx = 0.130;
INPUT.load = [   66  1 Fx/2
                 77  1 Fx/2
                 121 1 Fx
                 110 1 Fx
                 99  1 Fx
                 88  1 Fx ];

% -- Boundary conditions

spc_top = [ 61 2
62 2
63 2
64 2
65 2
66 2 ];

INPUT.spc = [   1 1
                7 1
               13 1
               19 1
               25 1
               31 1
               37 1
               43 1
               49 1
               55 1
               61 1
               67 2
               68 2
               69 2
               70 2
               71 2
               72 2
               73 2
               74 2
               75 2
               76 2
               77 2 
               spc_top ];

% --- Nonlinear solver parameters                                         
INPUT.sol_type     = 'nonlinear'; 
INPUT.lambda_max   = 1.0;
INPUT.dlambda      = 0.1;
INPUT.norm_res_max = 1e-10;
INPUT.nitermax     = 99;

