:MIXT LSO 7.4 3
  G4_Lu 0.764016
  G4_Si 0.0613283
  G4_O 0.174656

:ELEM Lutetium Lu 71 174.967

:MATE vacuum 1. 1. 1.e-25*g/cm3

// ROTATION MATRIX
:ROTM RM0 0. 0. 0. 

// THE WORLD
:VOLU world BOX 45.*mm 10.*mm 10.*mm G4_AIR

// SUB-WORLD
// box of air to put crystal blocks in
:VOLU CrystalAirBox BOX 10.*mm 6.*mm 6.*mm G4_AIR

// PLACEMENT OF SUB-WORLDS IN WORLD
// place 2 BOXES in world (front-faces) spaced 100 mm apart
// half crystal length (10 mm) plus distance from origin (20 mm)
:PLACE CrystalAirBox 1 world RM0  30.*mm 0.*mm 0.*mm
:PLACE CrystalAirBox 2 world RM0 -30.*mm 0.*mm 0.*mm

// PARAMETERS: CRYSTAL DIMENSIONS
// define dimensions of crystals in terms
// of half-lengths
:P crysYZ 2.*mm
:P crysX 10.*mm

// MOTHER 
// define volume as box which is a 3x3 array of crystals 
:VOLU CrystalBlock BOX $crysX 3*$crysYZ 3*$crysYZ LSO

// PLACEMENT OF MOTHER IN SUB-WORLD	
// place crystal blocks in AIR BOXES
:PLACE_PARAM CrystalBlock 1 CrystalAirBox SQUARE_YZ RM0 1 1 0.*mm 0.*mm 0.*mm 0.*mm

// DAUGHTER
// define volume as box which is a single crystal
:VOLU Crystal BOX  $crysX $crysYZ $crysYZ LSO

// PLACEMENT OF DAUGHTERS IN MOTHERS
// define how crystals will be placed in blocks 
:PLACE_PARAM Crystal 1 CrystalBlock SQUARE_YZ RM0 3 3 4*mm 4*mm -4*mm -4*mm

// USE AS EXTENDED SOURCE
// positron range plus source size
//:VOLU AirSphere ORB 2.*mm G4_AIR
//:PLACE AirSphere 1 world RM0 0. 0. 0. 
