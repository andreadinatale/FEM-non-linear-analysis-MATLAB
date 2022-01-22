function INPUT = inputNL_cantilever
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
INPUT.integration_pts = 3;
% INPUT.integration_barlow = 1; %optional field: to recover stresses with Barlow set to 1; otherwise omit or set to 0;

INPUT.elements = [  1      3     20     18      2     13     19     12
                    3      5     22     20      4     14     21     13
                    5      7     24     22      6     15     23     14
                    7      9     26     24      8     16     25     15
                    9     11     28     26     10     17     27     16  ];

% -- Nodes
INPUT.nodes = [   0.          0.  
  1.          0.  
  2.          0.  
  3.          0.  
  4.          0.  
  5.          0.  
  6.          0.  
  7.          0.  
  8.          0.  
  9.          0.  
 10.          0.  
  0.         0.5  
  2.         0.5  
  4.         0.5  
  6.         0.5  
  8.         0.5  
 10.         0.5  
  0.          1.  
  1.          1.  
  2.          1.  
  3.          1.  
  4.          1.  
  5.          1.  
  6.          1.  
  7.          1.  
  8.          1.  
  9.          1.  
 10.          1.    ];

% --- Material properties
INPUT.mat_type = 'plane_stress'; 
INPUT.E  = 1.2e7;
INPUT.nu = .3;
INPUT.t  = 0.1;
  
% -- Loading conditions
F0 = 1;
INPUT.load = [     1    2  -F0*2/12
                  11    2  -F0*2/12
                  28    2  -F0*2/12
                  18    2  -F0*2/12
                   2    2  -F0*2/3
                  19    2  -F0*2/3
                   4    2  -F0*2/3
                  21    2  -F0*2/3
                   6    2  -F0*2/3
                  23    2  -F0*2/3
                   8    2  -F0*2/3
                  25    2  -F0*2/3
                  10    2  -F0*2/3
                  27    2  -F0*2/3
                   3    2  -F0*2/6
                   5    2  -F0*2/6
                   7    2  -F0*2/6
                   9    2  -F0*2/6
                  26    2  -F0*2/6
                  24    2  -F0*2/6
                  22    2  -F0*2/6
                  20    2  -F0*2/6 ];

% -- Boundary conditions
INPUT.spc = [     1    1
    1    2
   12    1
   12    2
   18    1
   18    2 ];

% --- Nonlinear solver parameters                                         
INPUT.sol_type     = 'nonlinear'; 
INPUT.lambda_max   = 900.0;
INPUT.dlambda      = 50.;
INPUT.norm_res_max = 1e-10;
INPUT.nitermax     = 99;

