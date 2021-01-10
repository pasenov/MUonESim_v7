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
// $Id: EventAction.cc 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"
#include <iomanip>

#include <stdio.h>
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* DA, RunAction* RA)
:G4UserEventAction(), fDetectorconstruction(DA), fRunAction(RA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{

 fEnergyExitTargetPrimary = fEnergyExitTargetSecondaryPositron = fEnergyExitTargetSecondaryElectron = fEnergyExitTargetSecondaryMuonPlus = fEnergyExitTargetSecondaryMuonMinus = 0.;
 fAnglePrimary = fAngleSecondaryPositron = fAngleSecondaryElectron = fAngleSecondaryMuonPlus = fAngleSecondaryMuonMinus = 0.;

 fBoolPrim = fBoolSecPositron = fBoolSecElectron = fBoolSecMuonPlus = fBoolSecMuonMinus = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 // get event ID
 //G4int evtNb = evt->GetEventID();
 //G4cout 
    //<< "\n Event ID = " 
    //<< evtNb
    //<< G4endl;

 /*fRunAction->AddREnergyExitTargetPrimary(fEnergyExitTargetPrimary);
 fRunAction->AddREnergyExitTargetSecondaryPositron(fEnergyExitTargetSecondaryPositron);
 fRunAction->AddREnergyExitTargetSecondaryElectron(fEnergyExitTargetSecondaryElectron);
 fRunAction->AddREnergyExitTargetMuonPlus(fEnergyExitTargetMuonPlus);
 fRunAction->AddREnergyExitTargetMuonMinus(fEnergyExitTargetMuonMinus);

 fRunAction->AddRAnglePrimary(fAnglePrimary);
 fRunAction->AddRAngleSecondaryPositron(fAngleSecondaryPositron);
 fRunAction->AddRAngleSecondaryElectron(fAngleSecondaryElectron);
 fRunAction->AddRAngleMuonPlus(fAngleMuonPlus);
 fRunAction->AddRAngleMuonMinus(fAngleMuonMinus);*/

 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 analysisManager->SetVerboseLevel(1);
 
 if (fBoolPrim == 1)  {
    analysisManager->FillH1(1, fEnergyExitTargetPrimary);
    analysisManager->FillH1(6, fAnglePrimary);
 }
 if (fBoolSecPositron == 1)  {
    analysisManager->FillH1(2, fEnergyExitTargetSecondaryPositron);
    analysisManager->FillH1(7, fAngleSecondaryPositron);
 }
 if (fBoolSecElectron == 1)  {
    analysisManager->FillH1(3, fEnergyExitTargetSecondaryElectron);
    analysisManager->FillH1(8, fAngleSecondaryElectron);
 }
 if (fBoolSecMuonPlus == 1)  {
    analysisManager->FillH1(4, fEnergyExitTargetSecondaryMuonPlus);
    analysisManager->FillH1(9, fAngleSecondaryMuonPlus);
 }
 if (fBoolSecMuonMinus == 1)  {
    analysisManager->FillH1(5, fEnergyExitTargetSecondaryMuonMinus);
    analysisManager->FillH1(10, fAngleSecondaryMuonMinus);
 }

 if ((fBoolPrim == 1) && (fBoolSecPositron == 1))  {
    analysisManager->FillH2(1, fAnglePrimary, fAngleSecondaryPositron);
 }
 if ((fBoolPrim == 1) && (fBoolSecElectron == 1))  {
    analysisManager->FillH2(2, fAnglePrimary, fAngleSecondaryElectron);
    /*G4cout 
       << "\n Filling histogram 2."
       << G4endl;*/
 }
 if ((fBoolPrim == 1) && (fBoolSecMuonPlus == 1))  {
    analysisManager->FillH2(3, fAnglePrimary, fAngleSecondaryMuonPlus);
 }
 if ((fBoolPrim == 1) && (fBoolSecMuonMinus == 1))  {
    analysisManager->FillH2(4, fAnglePrimary, fAngleSecondaryMuonMinus);
 }

 //Visualize event if there is a condition
 //if (condition)  {
     //G4cout 
         //<< [condition] fullfilled
         //<< G4endl;
     //G4EventManager* evMan = G4EventManager::GetEventManager();
     //evMan->KeepTheCurrentEvent();
 //}

 //G4cout 
    //<< "\n End of event! "
    //<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

