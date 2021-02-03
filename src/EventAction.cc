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

 fEnergyExitTargetPrimary = fEnergyExitTargetSecondaryPositron = fEnergyExitTargetSecondaryElectron = fEnergyExitTargetSecondaryMuonPlus = fEnergyExitTargetSecondaryMuonMinus = fEnergyExitTargetCharged = fEnergyExitTargetChargedCut05GeV = fEnergyExitTargetChargedCut1GeV = fEnergyExitTargetIoni = fEnergyExitTargetPP1 = fEnergyExitTargetPP2 = fEnergyExitTargetIoniCut05GeV = fEnergyExitTargetPP1Cut05GeV = fEnergyExitTargetPP2Cut05GeV = fEnergyExitTargetIoniCut1GeV = fEnergyExitTargetPP1Cut1GeV = fEnergyExitTargetPP2Cut1GeV = 0.;
 fAnglePrimary = fAngleSecondaryPositron = fAngleSecondaryElectron = fAngleSecondaryMuonPlus = fAngleSecondaryMuonMinus = fAngleCharged = fAngleChargedCut05GeV = fAngleChargedCut1GeV = fAngleIoni = fAnglePP1 = fAnglePP2 = fAngleIoniCut05GeV = fAnglePP1Cut05GeV = fAnglePP2Cut05GeV = fAngleIoniCut1GeV = fAnglePP1Cut1GeV = fAnglePP2Cut1GeV = 0.;

 fBoolPrim = fBoolSecPositron = fBoolSecElectron = fBoolSecMuonPlus = fBoolSecMuonMinus = fBoolCharged = fBoolChargedCut05GeV = fBoolChargedCut1GeV = fBooleIoni = fBoolePP1 = fBoolePP2 = fBooleIoniCut05GeV = fBoolePP1Cut05GeV = fBoolePP2Cut05GeV = fBooleIoniCut1GeV = fBoolePP1Cut1GeV = fBoolePP2Cut1GeV = 0;
 fNbTracks = fNbTracksCut05GeV = fNbTracksCut1GeV = 0;

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
    analysisManager->FillNtupleDColumn(1, fEnergyExitTargetPrimary);
    analysisManager->FillNtupleDColumn(6, fAnglePrimary);
 }
 if (fBoolSecPositron == 1)  {
    analysisManager->FillH1(2, fEnergyExitTargetSecondaryPositron);
    analysisManager->FillH1(7, fAngleSecondaryPositron);
    analysisManager->FillNtupleDColumn(2, fEnergyExitTargetSecondaryPositron);
    analysisManager->FillNtupleDColumn(7, fAngleSecondaryPositron);
 }
 if (fBoolSecElectron == 1)  {
    analysisManager->FillH1(3, fEnergyExitTargetSecondaryElectron);
    analysisManager->FillH1(8, fAngleSecondaryElectron);
    analysisManager->FillNtupleDColumn(3, fEnergyExitTargetSecondaryElectron);
    analysisManager->FillNtupleDColumn(8, fAngleSecondaryElectron);
 }
 if (fBoolSecMuonPlus == 1)  {
    analysisManager->FillH1(4, fEnergyExitTargetSecondaryMuonPlus);
    analysisManager->FillH1(9, fAngleSecondaryMuonPlus);
    analysisManager->FillNtupleDColumn(4, fEnergyExitTargetSecondaryMuonPlus);
    analysisManager->FillNtupleDColumn(9, fAngleSecondaryMuonPlus);
 }
 if (fBoolSecMuonMinus == 1)  {
    analysisManager->FillH1(5, fEnergyExitTargetSecondaryMuonMinus);
    analysisManager->FillH1(10, fAngleSecondaryMuonMinus);
    analysisManager->FillNtupleDColumn(5, fEnergyExitTargetSecondaryMuonMinus);
    analysisManager->FillNtupleDColumn(10, fAngleSecondaryMuonMinus);
 }
 if (fBoolCharged == 1)  {
    analysisManager->FillH1(12, fEnergyExitTargetCharged);
    analysisManager->FillH1(13, fAngleCharged);
    analysisManager->FillNtupleDColumn(12, fEnergyExitTargetCharged);
    analysisManager->FillNtupleDColumn(13, fAngleCharged);
    
    analysisManager->FillH1(18, fNbTracks);
    analysisManager->FillNtupleDColumn(18, fNbTracks);
 }
 if (fBoolChargedCut05GeV == 1)  {
    analysisManager->FillH1(14, fEnergyExitTargetChargedCut05GeV);
    analysisManager->FillH1(15, fAngleChargedCut05GeV);
    analysisManager->FillNtupleDColumn(14, fEnergyExitTargetChargedCut05GeV);
    analysisManager->FillNtupleDColumn(15, fAngleChargedCut05GeV);
    
    analysisManager->FillH1(19, fNbTracksCut05GeV);
    analysisManager->FillNtupleDColumn(19, fNbTracksCut05GeV);
 }
 if (fBoolChargedCut1GeV == 1)  {
    analysisManager->FillH1(16, fEnergyExitTargetChargedCut1GeV);
    analysisManager->FillH1(17, fAngleChargedCut1GeV);
    analysisManager->FillNtupleDColumn(16, fEnergyExitTargetChargedCut1GeV);
    analysisManager->FillNtupleDColumn(17, fAngleChargedCut1GeV);
    
    analysisManager->FillH1(20, fNbTracksCut1GeV);
    analysisManager->FillNtupleDColumn(20, fNbTracksCut1GeV);
 }
 
 if (fBooleIoni == 1)  {
    analysisManager->FillH1(21, fEnergyExitTargetIoni);
    analysisManager->FillH1(27, fAngleIoni);
    
    analysisManager->FillNtupleDColumn(21, fEnergyExitTargetIoni);
    analysisManager->FillNtupleDColumn(27, fAngleIoni);
 }
 if (fBoolePP1 == 1)  {
    analysisManager->FillH1(22, fEnergyExitTargetPP1);
    analysisManager->FillH1(28, fAnglePP1);
    
    analysisManager->FillNtupleDColumn(22, fEnergyExitTargetPP1);
    analysisManager->FillNtupleDColumn(28, fAnglePP1);
 }
 if (fBoolePP2 == 1)  {
    analysisManager->FillH1(22, fEnergyExitTargetPP2);
    analysisManager->FillH1(28, fAnglePP2);
    
    analysisManager->FillNtupleDColumn(22, fEnergyExitTargetPP2);
    analysisManager->FillNtupleDColumn(28, fAnglePP2);
 }
 
 if (fBooleIoniCut05GeV == 1)  {
    analysisManager->FillH1(23, fEnergyExitTargetIoniCut05GeV);
    analysisManager->FillH1(29, fAngleIoniCut05GeV);
    
    analysisManager->FillNtupleDColumn(23, fEnergyExitTargetIoniCut05GeV);
    analysisManager->FillNtupleDColumn(29, fAngleIoniCut05GeV);
 }
 if (fBoolePP1Cut05GeV == 1)  {
    analysisManager->FillH1(24, fEnergyExitTargetPP1Cut05GeV);
    analysisManager->FillH1(30, fAnglePP1Cut05GeV);
    
    analysisManager->FillNtupleDColumn(24, fEnergyExitTargetPP1Cut05GeV);
    analysisManager->FillNtupleDColumn(30, fAnglePP1Cut05GeV);
 }
 if (fBoolePP2Cut05GeV == 1)  {
    analysisManager->FillH1(24, fEnergyExitTargetPP2Cut05GeV);
    analysisManager->FillH1(30, fAnglePP2Cut05GeV);
    
    analysisManager->FillNtupleDColumn(24, fEnergyExitTargetPP2Cut05GeV);
    analysisManager->FillNtupleDColumn(30, fAnglePP2Cut05GeV);
 }
 
 if (fBooleIoniCut1GeV == 1)  {
    analysisManager->FillH1(25, fEnergyExitTargetIoniCut1GeV);
    analysisManager->FillH1(31, fAngleIoniCut1GeV);
    
    analysisManager->FillNtupleDColumn(25, fEnergyExitTargetIoniCut1GeV);
    analysisManager->FillNtupleDColumn(31, fAngleIoniCut1GeV);
 }
 if (fBoolePP1Cut1GeV == 1)  {
    analysisManager->FillH1(26, fEnergyExitTargetPP1Cut1GeV);
    analysisManager->FillH1(32, fAnglePP1Cut1GeV);
    
    analysisManager->FillNtupleDColumn(26, fEnergyExitTargetPP1Cut1GeV);
    analysisManager->FillNtupleDColumn(32, fAnglePP1Cut1GeV);
 }
 if (fBoolePP2Cut1GeV == 1)  {
    analysisManager->FillH1(26, fEnergyExitTargetPP2Cut1GeV);
    analysisManager->FillH1(32, fAnglePP2Cut1GeV);
    
    analysisManager->FillNtupleDColumn(26, fEnergyExitTargetPP2Cut1GeV);
    analysisManager->FillNtupleDColumn(32, fAnglePP2Cut1GeV);
 }

 if ((fBoolPrim == 1) && (fBoolSecPositron == 1))  {
    analysisManager->FillH2(1, fAngleSecondaryPositron, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolSecElectron == 1))  {
    analysisManager->FillH2(2, fAngleSecondaryElectron, fAnglePrimary);
    /*G4cout 
       << "\n Filling histogram 2."
       << G4endl;*/
 }
 if ((fBoolPrim == 1) && (fBoolSecMuonPlus == 1))  {
    analysisManager->FillH2(3, fAngleSecondaryMuonPlus, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolSecMuonMinus == 1))  {
    analysisManager->FillH2(4, fAngleSecondaryMuonMinus, fAnglePrimary);
 }
 
 if ((fBoolPrim == 1) && (fBooleIoni == 1))  {
    analysisManager->FillH2(5, fAngleIoni, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolePP1 == 1))  {
    analysisManager->FillH2(6, fAnglePP1, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolePP2 == 1))  {
    analysisManager->FillH2(6, fAnglePP2, fAnglePrimary);
 }
 
 if ((fBoolPrim == 1) && (fBooleIoniCut05GeV == 1))  {
    analysisManager->FillH2(7, fAngleIoniCut05GeV, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolePP1Cut05GeV == 1))  {
    analysisManager->FillH2(8, fAnglePP1Cut05GeV, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolePP2Cut05GeV == 1))  {
    analysisManager->FillH2(8, fAnglePP2Cut05GeV, fAnglePrimary);
 }
 
 if ((fBoolPrim == 1) && (fBooleIoniCut1GeV == 1))  {
    analysisManager->FillH2(9, fAngleIoniCut1GeV, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolePP1Cut1GeV == 1))  {
    analysisManager->FillH2(10, fAnglePP1Cut1GeV, fAnglePrimary);
 }
 if ((fBoolPrim == 1) && (fBoolePP2Cut1GeV == 1))  {
    analysisManager->FillH2(10, fAnglePP2Cut1GeV, fAnglePrimary);
 }
 
 
 
 analysisManager->AddNtupleRow();

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

