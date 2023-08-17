La = 0.001;
L = 0.01;
Lb = 0.001;
h = 0.00125;
lc = L/60.; //0.1*h;
delta = 1e-8; //boundary layer 1
delta2 = 1e-6;// refined region around points
delta2b = 2e-4; //start of electrode
delta3 = 0.4*h; //boundary layer 2
delta3b = 0.01*h; //transition between layers
L_defect = .001; //length of defect
lc2 = lc;

Point(1) = {0, 0, 0, lc};
Point(2) = {La-delta2, 0, 0, lc2};
Point(3) = {La, 0, 0, lc2};
Point(4) = {La+0.1*delta2b, 0, 0, lc2};
Point(5) = {La+delta2b, 0, 0, lc2};
Point(6) = {La+L/2-delta2, 0, 0, lc2};
Point(7) = {La+L/2, 0, 0, lc2};
Point(8) = {La+L/2+delta2, 0, 0, lc2};
Point(9) = {La+L/2+L_defect-delta2, 0, 0, lc2};
Point(10) = {La+L/2+L_defect, 0, 0, lc2};
Point(11) = {La+L/2+L_defect+0.1*delta2b, 0, 0, lc2};
Point(12) = {La+L/2+L_defect+delta2b, 0, 0, lc2};
Point(13) = {La+L/2+L_defect+L/2-delta2, 0, 0, lc2};
Point(14) = {La+L/2+L_defect+L/2, 0, 0, lc2};
Point(15) = {La+L/2+L_defect+L/2+delta2, 0, 0, lc2};
Point(16) = {La+L/2+L_defect+L/2+Lb, 0, 0, lc};

Point(17) = {0, delta, 0, lc2};
Point(18) = {La-delta2, delta, 0, lc2};
Point(19) = {La, delta, 0, lc2};
Point(20) = {La+0.1*delta2b, delta, 0, lc2};
Point(21) = {La+delta2b, delta, 0, lc2};
Point(22) = {La+L/2-delta2, delta, 0, lc2};
Point(23) = {La+L/2, delta, 0, lc2};
Point(24) = {La+L/2+delta2, delta, 0, lc2};
Point(25) = {La+L/2+L_defect-delta2, delta, 0, lc2};
Point(26) = {La+L/2+L_defect, delta, 0, lc2};
Point(27) = {La+L/2+L_defect+0.1*delta2b, delta, 0, lc2};
Point(28) = {La+L/2+L_defect+delta2b, delta, 0, lc2};
Point(29) = {La+L/2+L_defect+L/2-delta2, delta, 0, lc2};
Point(30) = {La+L/2+L_defect+L/2, delta, 0, lc2};
Point(31) = {La+L/2+L_defect+L/2+delta2, delta, 0, lc2};
Point(32) = {La+L/2+L_defect+L/2+Lb, delta, 0, lc2};


Point(33) = {0, h, 0, lc};
Point(34) = {La+L/2+L/2+L_defect+Lb, h, 0, lc};

Point(35) = {La + delta2b, delta3b, 0, lc2};
Point(36) = {La + L/2 - delta2, delta3b, 0, lc2};

Point(37) = {La + delta2b, delta3, 0, lc2};
Point(38) = {La + L/2 - delta2, delta3, 0, lc2};

Point(39) = {La + L/2 + L_defect + delta2b, delta3b, 0, lc2};
Point(40) = {La + L/2 + L_defect + L/2 - delta2, delta3b, 0, lc2};

Point(41) = {La + L/2 + L_defect + delta2b, delta3, 0, lc2};
Point(42) = {La + L/2 + L_defect + L/2 - delta2, delta3, 0, lc2};


//horizontal lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};

Line(16) = {17,18};
Line(17) = {18,19};
Line(18) = {19,20};
Line(19) = {20,21};
Line(20) = {21,22};
Line(21) = {22,23};
Line(22) = {23,24};
Line(23) = {24,25};
Line(24) = {25,26};
Line(25) = {26,27};
Line(26) = {27,28};
Line(27) = {28,29};
Line(28) = {29,30};
Line(29) = {30,31};
Line(30) = {31,32};

Line(31) = {33,34};

//vertical lines
Line(32) = {1,17};
Line(33) = {2,18};
Line(34) = {3,19};
Line(35) = {4,20};
Line(36) = {5,21};
Line(37) = {6,22};
Line(38) = {7,23};
Line(39) = {8,24};
Line(40) = {9,25};
Line(41) = {10,26};
Line(42) = {11,27};
Line(43) = {12,28};
Line(44) = {13,29};
Line(45) = {14,30};
Line(46) = {15,31};
Line(47) = {16,32};

Line(48) = {17,33};
Line(49) = {32,34};

// extra lines
Line(50) = {35,36};
Line(51) = {21,35};
Line(52) = {22,36};

Line(53) = {37,38};
Line(54) = {35,37};
Line(55) = {36,38};

Line(56) = {39,40};
Line(57) = {28,39};
Line(58) = {29,40};

Line(59) = {41,42};
Line(60) = {39,41};
Line(61) = {40,42};

Curve Loop(1) = {1,33,-16,-32};
Curve Loop(2) = {2,34,-17,-33};
Curve Loop(3) = {3,35,-18,-34};
Curve Loop(4) = {4,36,-19,-35};
Curve Loop(5) = {5,37,-20,-36};
Curve Loop(6) = {6,38,-21,-37};
Curve Loop(7) = {7,39,-22,-38};
Curve Loop(8) = {8,40,-23,-39};
Curve Loop(9) = {9,41,-24,-40};
Curve Loop(10) = {10,42,-25,-41};
Curve Loop(11) = {11,43,-26,-42};
Curve Loop(12) = {12,44,-27,-43};
Curve Loop(13) = {13,45,-28,-44};
Curve Loop(14) = {14,46,-29,-45};
Curve Loop(15) = {15,47,-30,-46};

Curve Loop(16) = {16,17,18,19,51,54,53,-55,-52,21,22,23,24,25,26,57,60,59,-61,-58,28,29,30,49,-31,-48};

Curve Loop(17) = {20,52,-50,-51};
Curve Loop(18) = {50,55,-53,-54};

Curve Loop(19) = {27,58,-56,-57};
Curve Loop(20) = {56,61,-59,-60};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};

Plane Surface(16) = {17};
Plane Surface(17) = {18};
Plane Surface(18) = {19};
Plane Surface(19) = {20};

//the unstructred part is meshed last
Plane Surface(20) = {16};

Physical Curve("Inlet", 1) = {32,48};
Physical Curve("Outlet", 2) = {47,49};
Physical Curve("Bulk", 3) = {31};
Physical Curve("Cathode", 4) = {3,4,5,6,10,11,12,13};
Physical Curve("Wall", 5) = {1,2,7,8,9,14,15};
Physical Surface("Channel", 1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

//y-dir
Transfinite Curve{32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47} = 26 Using Progression 1.15;

//x-dir electrode
Transfinite Curve{2,17,9,24} = 29 Using Progression 1/1.4;
Transfinite Curve{3,18,10,25} = 41 Using Progression 1.35;
Transfinite Curve{4,19,11,26} = 25 Using Progression 1.0;
Transfinite Curve{6,21,13,28} = 29 Using Progression 1/1.4;
Transfinite Curve{7,22,14,29} = 29 Using Progression 1.4;
Transfinite Curve{5,20,50,53,12,27,56,59} = 49 Using Bump 0.4;

//x-dir inlet transition
Transfinite Curve{1,16} = 51 Using Progression .9;

//x-dir outlet transition
Transfinite Curve{15,30} = 71 Using Progression 1./.95;

//x-dir defect
Transfinite Curve{8,23} = 31 Using Bump 0.1;

//y-dir diffusion bdlayer
Transfinite Curve{51,52,57,58} = 7 Using Progression 1.4;
Transfinite Curve{54,55,60,61} = 21 Using Progression 1;

//boundary layer
Transfinite Surface{1} = {1,2,18,17};
Transfinite Surface{2} = {2,3,19,18};
Transfinite Surface{3} = {3,4,20,19};
Transfinite Surface{4} = {4,5,21,20};
Transfinite Surface{5} = {5,6,22,21};
Transfinite Surface{6} = {6,7,23,22};
Transfinite Surface{7} = {7,8,24,23};
Transfinite Surface{8} = {8,9,25,24};
Transfinite Surface{9} = {9,10,26,25};
Transfinite Surface{10} = {10,11,27,26};
Transfinite Surface{11} = {11,12,28,27};
Transfinite Surface{12} = {12,13,29,28};
Transfinite Surface{13} = {13,14,30,29};
Transfinite Surface{14} = {14,15,31,30};
Transfinite Surface{15} = {15,16,32,31};

Transfinite Surface{16} = {21,22,36,35};
Transfinite Surface{17} = {35,36,38,37};
Transfinite Surface{18} = {28,29,40,39};
Transfinite Surface{19} = {39,40,42,41};


// Electrode start and end
Field[1] = Distance;
// Field[1].PointsList = {10,11,12,13,14,15};
Field[1].PointsList = {19,23,26,30};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 8e-11;
Field[2].SizeMax = lc;
Field[2].DistMin = 2e-10;
Field[2].DistMax = 0.8*h;

Field[3] = Distance;
Field[3].CurvesList = {20,21,22,23,27,28,29,30};
Field[3].NumPointsPerCurve = 100;

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 0.01/100.;
Field[4].SizeMax = lc;
Field[4].DistMin = 0.5*h;
Field[4].DistMax = h;
//Field[4].Sigmoid = 1;

Field[5] = Distance;
Field[5].CurvesList = {30};
Field[5].NumPointsPerCurve = 100;

Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = 0.3*lc;
Field[6].SizeMax = lc;
Field[6].DistMin = delta3;
Field[6].DistMax = 0.5*h;

Field[7] = Distance;
Field[7].CurvesList = {17,24};
Field[7].NumPointsPerCurve = 100;

Field[8] = Threshold;
Field[8].InField = 7;
Field[8].SizeMin = 0.8*delta2b/22;
Field[8].SizeMax = lc;
Field[8].DistMin = delta3*0.1;
Field[8].DistMax = 0.5*h;

Field[9] = Min;
Field[9].FieldsList = {2,4,6, 8}; 
Background Field = 9;

Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
//Recombine Surface{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
Mesh.RecombineAll=1;
