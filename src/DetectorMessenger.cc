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
/// \file electromagnetic/TestEm18/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),    
 fMaterWCmd(0),
 fMaterTCmd(0),
 fSizeXWCmd(0),
 fSizeYWCmd(0),
 fSizeZWCmd(0),
 fSizeXTCmd(0),
 fSizeYTCmd(0),
 fSizeZTCmd(0),
 fUpdateCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance("commands specific to this example");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction");
        
  fMaterWCmd = new G4UIcmdWithAString("/testem/det/setMatW",this);
  fMaterWCmd->SetGuidance("Select material of the world.");
  fMaterWCmd->SetParameterName("choice",false);
  fMaterWCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMaterTCmd = new G4UIcmdWithAString("/testem/det/setMatT",this);
  fMaterTCmd->SetGuidance("Select material of the target.");
  fMaterTCmd->SetParameterName("choice",false);
  fMaterTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeXWCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeXW",this);
  fSizeXWCmd->SetGuidance("Set size-X of the world");
  fSizeXWCmd->SetParameterName("SizeXW",false);
  fSizeXWCmd->SetRange("SizeXW>0.");
  fSizeXWCmd->SetUnitCategory("Length");
  fSizeXWCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeYWCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeYW",this);
  fSizeYWCmd->SetGuidance("Set size-Y of the world");
  fSizeYWCmd->SetParameterName("SizeYW",false);
  fSizeYWCmd->SetRange("SizeYW>0.");
  fSizeYWCmd->SetUnitCategory("Length");
  fSizeYWCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSizeZWCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeZW",this);
  fSizeZWCmd->SetGuidance("Set size-Z of the world");
  fSizeZWCmd->SetParameterName("SizeZW",false);
  fSizeZWCmd->SetRange("SizeZW>0.");
  fSizeZWCmd->SetUnitCategory("Length");
  fSizeZWCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeXTCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeXT",this);
  fSizeXTCmd->SetGuidance("Set size-X of the target");
  fSizeXTCmd->SetParameterName("SizeXT",false);
  fSizeXTCmd->SetRange("SizeXT>0.");
  fSizeXTCmd->SetUnitCategory("Length");
  fSizeXTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeYTCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeYT",this);
  fSizeYTCmd->SetGuidance("Set size-Y of the target");
  fSizeYTCmd->SetParameterName("SizeYT",false);
  fSizeYTCmd->SetRange("SizeYT>0.");
  fSizeYTCmd->SetUnitCategory("Length");
  fSizeYTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeZTCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeZT",this);
  fSizeZTCmd->SetGuidance("Set size-Z of the target");
  fSizeZTCmd->SetParameterName("SizeZT",false);
  fSizeZTCmd->SetRange("SizeZT>0.");
  fSizeZTCmd->SetUnitCategory("Length");
  fSizeZTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  fUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterWCmd;
  delete fMaterTCmd;
  delete fSizeXWCmd;
  delete fSizeYWCmd;
  delete fSizeZWCmd;
  delete fSizeXTCmd; 
  delete fSizeYTCmd; 
  delete fSizeZTCmd; 
  delete fUpdateCmd;
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fMaterWCmd )
   { fDetector->SetMaterialW(newValue);}

  if( command == fMaterTCmd )
   { fDetector->SetMaterialT(newValue);}

  if( command == fSizeXWCmd )
   { fDetector->SetSizeXW(fSizeXWCmd->GetNewDoubleValue(newValue));}

  if( command == fSizeYWCmd )
   { fDetector->SetSizeYW(fSizeYWCmd->GetNewDoubleValue(newValue));}
   
  if( command == fSizeZWCmd )
   { fDetector->SetSizeZW(fSizeZWCmd->GetNewDoubleValue(newValue));}

  if( command == fSizeXTCmd )
   { fDetector->SetSizeXT(fSizeXTCmd->GetNewDoubleValue(newValue));}

  if( command == fSizeYTCmd )
   { fDetector->SetSizeYT(fSizeYTCmd->GetNewDoubleValue(newValue));}

  if( command == fSizeZTCmd )
   { fDetector->SetSizeZT(fSizeZTCmd->GetNewDoubleValue(newValue));}
     
  if( command == fUpdateCmd )
   { fDetector->UpdateGeometry(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
