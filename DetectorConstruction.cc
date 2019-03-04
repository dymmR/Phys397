//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr03/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double density;
G4double a;
DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fTargetMater(0), fLogicTarget(0),
 fDetectorMater(0), fLogicDetector(0),
 fWorldMater(0), fPhysiWorld(0), fCupMater(0),
 fLidMaterial(0),
 fDetectorMessenger(0)
{
  //fTargetLength      = 1*cm;
  //fTargetRadius      = 0.5*cm;
  fDetectorLength    = 4.31*cm;
  fDetectorThickness = 2.77*cm;

  fWorldLength = std::max(fTargetLength,fDetectorLength);
  fWorldRadius = fTargetRadius + fDetectorThickness;

  DefineMaterials();

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // build materials
  //
  fDetectorMater =
  new G4Material("Germanium", 32, 72.61*g/mole, 5.323*g/cm3);


  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
  //

  G4int ncomponents; G4double fractionmass;
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
                      kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
    fTargetMater=new G4Material("Copper", 29,63.546*g/mole, 8.96*g/cm3);
  //
  fWorldMater = Air20;



  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  fCupMater = new G4Material("Aluminum", 13., a, density);
    fLidMaterial = new G4Material("Lead", 82., 207.2*g/mole, 11.34*g/cm3);

  // or use G4 materials data base
  //
  G4NistManager* man = G4NistManager::Instance();


 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // World
  //
  // (re) compute World dimensions if necessary


  G4double target_x=0.5*cm;
  G4double target_y=10*cm;
  G4double target_z=10*cm;

  fWorldLength = 25*cm;
  fWorldRadius = 25*cm;

  G4Tubs*
  sWorld = new G4Tubs("World",                                 //name
                 0.,fWorldRadius, fWorldLength, 0.,twopi); //dimensions

  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMater,               //material
                             "World");                  //name

  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                            G4ThreeVector(0,0,0),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number

  // Target
  //

  G4RotationMatrix* yRot = new G4RotationMatrix; // Rotates X and Z axes only
yRot-> rotateY(M_PI/4.*rad);

  G4Tubs* sTarget = new G4Tubs("target",0, 5*cm,  0.5*cm, 0., twopi); //dimensions

G4String target_name="target";
G4LogicalVolume*
fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "target");

                                           //name
G4VPhysicalVolume* fTargetPhys=new G4PVPlacement(yRot,                         //no rotation
                           G4ThreeVector(1.76*cm,0,27*cm),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           target_name,                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number


//old
//G4Tubs*
//sTarget = new G4Tubs("Target",                                   //name
  //              0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions


//fLogicTarget = new G4LogicalVolume(sTarget,           //shape
  //                         fTargetMater,              //material
    //                       "Target");                 //name

      //   new G4PVPlacement(rotateY(M_PI/4.*rad),                         //no rotation
        //                 G4ThreeVector(0,0,14.*cm),             //at (0,0,0)
          //               fLogicTarget,                //logical volume
            //             "Target",                    //name
              //           lWorld,                      //mother  volume
                //         false,                       //no boolean operation
                  //       0);


                G4Tubs*  sContainer = new G4Tubs("Container",                                 //name
                                 0,15*cm, 16*cm, 0,twopi); //dimensions
                                 //there was a 0.5*16

G4LogicalVolume*
                  fLogicContainer = new G4LogicalVolume(sContainer,                  //shape
                                             fLidMaterial,               //material
                                             "Container");


                                             G4String container_name="Container";             //name
G4VPhysicalVolume* fContainerPhys= new G4PVPlacement(0,                    //no rotation
                                            G4ThreeVector(),
                                            fLogicContainer,           //at (0,0,0)
                                                               //logical volume
                                            container_name,
                                            lWorld,                     //name
                                             //mother volume
                                            false,                      //no boolean operation
                                            0);





  // Detector
  //
  G4Tubs*
  sDetector = new G4Tubs("Detector",
                fTargetRadius, fWorldRadius, 0.5*fDetectorLength, 0.,twopi);

G4LogicalVolume*
  fLogicDetector = new G4LogicalVolume(sDetector,       //shape
                             fDetectorMater,            //material
                             "Detector");               //name
G4VPhysicalVolume* fDetectorPhys=
           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(0,0,0),             //at (0,0,0)
                           fLogicDetector,              //logical volume
                           "Detector",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number





//aluminum plating

 G4Tubs* side_aluminum = new G4Tubs("side_aluminum",
                    	3.4*cm,
                   		3.5*cm,
                    	2.455*cm,
                   		0,
                   		2*M_PI);

G4Tubs* top_aluminum = new G4Tubs("top_aluminum",
                              0*cm,
                           		3.5*cm,
                           		0.5*mm,
                           		0,
                           		2*M_PI);


G4UnionSolid* fullCup = new G4UnionSolid("Cup", side_aluminum, top_aluminum, NULL, G4ThreeVector(0,0,2.505*cm));

G4LogicalVolume*
         fLogicCup = new G4LogicalVolume(fullCup,       //shape
                                          fCupMater,            //material
                                          "Cup");               //name
G4VPhysicalVolume* fCupPhys=
            new G4PVPlacement(0,                         //no rotation
                                    G4ThreeVector(0,0,0),             //at (0,0,0)
                                      fLogicCup,              //logical volume
                                          "Cup",                  //name
                                          lWorld,                      //mother  volume
                                          false,                       //no boolean operation
                                          0);                          //copy number


//lead top

G4VSolid* lid=new G4Box("lid", 7*cm, 7*cm, 1.5*cm);
G4Tubs* inner_hole=new G4Tubs("lid hole",0, 0.35*cm,0.5*cm,0,2*M_PI);

G4SubtractionSolid *full_lid=new G4SubtractionSolid("full_lid", lid, inner_hole);

G4LogicalVolume*
fLogicLid=new G4LogicalVolume(full_lid, fLidMaterial, "full_lid");


G4VPhysicalVolume* fLidPhys=
new G4PVPlacement(0,                         //no rotation
                        G4ThreeVector(0,0,10.5*cm),             //at (0,0,0)
                          fLogicLid,              //logical volume
                              "full_lid",                  //name
                              lWorld,                      //mother  volume
                              false,                       //no boolean operation
                              0);                          //copy number





  G4VSolid*
          sBlocker = new G4Box("Blocker",                                   //name
                                              5.5*cm, 10*cm, 10.5*cm); //dimensions

G4LogicalVolume*
                              fLogicBlocker = new G4LogicalVolume(sBlocker,           //shape
                                                         fLidMaterial,              //material
                                                         "Blocker");                 //name
G4VPhysicalVolume* BlockerPhys=
                                       new G4PVPlacement(0,                         //no rotation
                                                       G4ThreeVector(5*cm,0,0),             //at (0,0,0)
                                                       fLogicBlocker,                //logical volume
                                                       "Blocker",                    //name
                                                       lWorld,                      //mother  volume
                                                       false,                       //no boolean operation
                                                       0);








  PrintParameters();

  //always return the root volume
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")
         << " Material = " << fTargetMater->GetName();
  G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
         << " Tickness = " << G4BestUnit(fDetectorThickness,"Length")
         << " Material = " << fDetectorMater->GetName() << G4endl;
  G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMater = pttoMaterial;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fDetectorMater = pttoMaterial;
    if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}



void DetectorConstruction::SetCupMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fCupMater = pttoMaterial;
    if(fLogicCup) { fLogicCup->SetMaterial(fCupMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}


void DetectorConstruction::SetLidMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fLidMaterial = pttoMaterial;
    if(fLogicLid) { fLogicLid->SetMaterial(fLidMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
  fTargetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
  fTargetLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
  fDetectorThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
  fDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
  return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
  return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
  return fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
  return fLogicDetector;
}




G4Material* DetectorConstruction::GetCupMaterial()
{
  return fCupMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicCup()
{
  return fLogicCup;
}

G4Material* DetectorConstruction::GetBlockerMaterial()
{
  return fLidMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicBlocker()
{
  return fLogicBlocker;
}

G4Material* DetectorConstruction::GetLidMaterial()
{
  return fLidMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicLid()
{
  return fLogicLid;
}


G4Material* DetectorConstruction::GetContainerMaterial()
{
  return fLidMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicContainer()
{
  return fLogicContainer;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
