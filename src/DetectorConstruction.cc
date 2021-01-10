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
// $Id: DetectorConstruction.cc 68348 2013-03-22 10:00:19Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ElectricFieldSetup.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <iomanip>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4AutoDelete.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4NistMaterialBuilder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fPBoxW(0), fPBoxT(0), fLBoxW(0), fLBoxT(0), fMaterialW(0), fMaterialT(0), fDetectorMessenger(0)
{
  zOff = 13.5*cm;

  fWorldSizeX   = 10.0*cm; 
  fWorldSizeY   = 10.0*cm; 
  fWorldSizeZ   = 10.0*cm;

  fTargetSizeX = 1.5*cm;
  fTargetSizeY = 1.5*cm;
  fTargetSizeZ = 1.5*cm;

  DefineMaterials();
  SetMaterialW("Air");
  SetMaterialT("Graphite"); //Or beryllium

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
  //
  // define Elements
  //
  G4double z,a;
  
  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"   ,"C" , z= 6., a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
  
  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;
  G4double fractionmass;  

  G4Material* H2O = 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  G4Material* CO2 = 
  new G4Material("CarbonicGas", density= 1.87*mg/cm3, ncomponents=2, kStateGas, 273.15*kelvin, 1*atmosphere);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);
  CO2->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  G4Material* vapor = 
  new G4Material("Water_vapor", density= 1.000*mg/cm3, ncomponents=2);
  vapor->AddElement(H, natoms=2);
  vapor->AddElement(O, natoms=1);
  vapor->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  new G4Material("Carbon"     , z=6.,  a= 12.01*g/mole, density= 2.267*g/cm3);
  new G4Material("Graphite"   , z=6.,  a= 12.00*g/mole, density= 2.265*g/cm3);
  //Graphite->SetChemicalFormula("Graphite");
  new G4Material("Aluminium"  , z=13., a= 26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"    , z=14., a= 28.09*g/mole, density= 2.330*g/cm3);
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  new G4Material("Iron"       , z=26., a= 55.85*g/mole, density= 7.870*g/cm3);  
  new G4Material("Germanium"  , z=32., a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Tungsten"   , z=74., a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Lead"       , z=82., a=207.19*g/mole, density= 11.35*g/cm3);
  new G4Material("Nitrogen"   , z=7.,  a= 14.01*g/mole, density= 1.145*mg/cm3);
  
  G4Material* ArgonGas =   
  new G4Material("ArgonGas"   , z=18., a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);


  G4Material* Air = 
  new G4Material("Air", density= 1.290*mg/cm3, ncomponents=4);
  Air->AddElement(N, fractionmass=78.08*perCent);
  Air->AddElement(O, fractionmass=20.95*perCent);
  Air->AddMaterial(ArgonGas, fractionmass=0.93*perCent);
  Air->AddMaterial(CO2, fractionmass=0.04*perCent);

  G4Material* Scintillator = 
  new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Scintillator->AddElement(C, natoms=9);
  Scintillator->AddElement(H, natoms=10);
  Scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4RotationMatrix* myRotationTarget = new G4RotationMatrix();
  myRotationTarget->rotateX(0.*deg);
  myRotationTarget->rotateY(0.*deg);
  myRotationTarget->rotateZ(0.*deg);

  G4Box*
  sBoxW = new G4Box("World",                              //its name
                   fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2);       //its dimensions

  fLBoxW = new G4LogicalVolume(sBoxW,                        //its shape
                             fMaterialW,                    //its material
                             fMaterialW->GetName());        //its name

  fPBoxW = new G4PVPlacement(0,                          //no rotation
                           G4ThreeVector(),           //at (0,0,0)
                           fLBoxW,                       //its logical volume
                           fMaterialW->GetName(),        //its name
                           0,                           //its mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // Target

  G4Box*
  sBoxT = new G4Box("Target",				     //its name
				 fTargetSizeX/2,fTargetSizeY/2,fTargetSizeZ/2);   //its dimensions
  

  fLBoxT = new G4LogicalVolume(sBoxT,                   //its shape
			     fMaterialT,                   //its material
			     fMaterialT->GetName());       //its name
  
  pos_xT =  0.0*cm;
  pos_yT =  0.0*cm;
  pos_zT =  0.0*cm;
  
  fPBoxT =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xT, pos_yT, pos_zT),	//at (pos_xT, pos_yT, pos_zT)
                    "Target",                  //its name
                    fLBoxT,                       //its logical volume
                    fPBoxW,                           //its mother  volume
                    false,                       //no boolean operation
                    0);                          //copy number

 

  // Visualization attributes
  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  worldVisAtt->SetVisibility(true);
  fLBoxW->SetVisAttributes(worldVisAtt);

  G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5)); //Grey
  targetVisAtt->SetVisibility(true);
  fLBoxT->SetVisAttributes(targetVisAtt);
                           
  PrintParameters();
  
  //always return the root volume
  //
  return fPBoxW;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The World is " << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fMaterialW->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialW(const G4String& nameW)
{
  // search the material by its name
  G4Material* matW = G4Material::GetMaterial(nameW, false);

  // create the material by its name
  if(!matW) { matW = G4NistManager::Instance()->FindOrBuildMaterial(nameW); }

  if(matW != fMaterialW) {
    G4cout << "### New material " << matW->GetName() << G4endl;
    fMaterialW = matW;
    UpdateGeometry();
  }

  if(!matW) {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterialW : "
           << nameW << " not found" << G4endl;  
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialT(const G4String& nameT)
{
  // search the material by its name
  G4Material* matT = G4Material::GetMaterial(nameT, false);

  // create the material by its name
  if(!matT) { matT = G4NistManager::Instance()->FindOrBuildMaterial(nameT); }

  if(matT != fMaterialT) {
    G4cout << "### New material " << matT->GetName() << G4endl;
    fMaterialT = matT;
    UpdateGeometry();
  }

  if(!matT) {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterialT : "
           << nameT << " not found" << G4endl;  
  } 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeXW(G4double valueXW)
{
  fWorldSizeX = valueXW;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeYW(G4double valueYW)
{
  fWorldSizeY = valueYW;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeZW(G4double valueZW)
{
  fWorldSizeZ = valueZW;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeXT(G4double valueXT)
{
  fTargetSizeX = valueXT;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeYT(G4double valueYT)
{
  fTargetSizeY = valueYT;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeZT(G4double valueZT)
{
  fTargetSizeZ = valueZT;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive Detectors Absorber

  //if (!fCalorimeterSD.Get()) {
    //CalorimeterSD* calorimeterSD = new CalorimeterSD("CalorSD",this);
    //fCalorimeterSD.Put(calorimeterSD);
  //}  
  //G4SDManager::GetSDMpointer()->AddNewDetector(fCalorimeterSD.Get());
  //SetSensitiveDetector(fLogicAbsorber, fCalorimeterSD.Get());

  // Construct the field creator - this will register the field it creates

  if (!fEmFieldSetup.Get()) { 
    ElectricFieldSetup* fieldSetup = new ElectricFieldSetup();
    G4AutoDelete::Register(fieldSetup); //Kernel will delete the messenger
    fEmFieldSetup.Put(fieldSetup);
  } 
 

  fElFieldz = 0;
  //G4cout 
      //<< "\n Electric field inside Target: " 
      //<< G4BestUnit(fElFieldz, "Electric field")
      //<< G4endl;

  fElField = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, 0));

  fLocalEquation = new G4EqMagElectricField(fElField);

  G4int nvar = 8;
  fLocalStepper = new G4ClassicalRK4(fLocalEquation, nvar);

  G4double fMinStep = 0.010*mm;

  fIntgrDriver = new G4MagInt_Driver(fMinStep, fLocalStepper, fLocalStepper->GetNumberOfVariables());
  fLocalChordFinder = new G4ChordFinder(fIntgrDriver);

  fLocalFieldManager = new G4FieldManager();
  fLocalFieldManager->SetDetectorField(fElField);
  fLocalFieldManager->SetChordFinder(fLocalChordFinder);

  G4bool allLocal = true ;
  fLBoxT->SetFieldManager(fLocalFieldManager, allLocal);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
