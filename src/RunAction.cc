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
/// \file electromagnetic/TestEm18/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin), fHistoManager(0)
{ 
  fHistoManager = new HistoManager(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
  delete fHistoManager; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //initialization
  //
  fEnergyExitTargetPrimary = 0.;
  fEnergyExitTargetSecondaryPositron = 0.;
  fEnergyExitTargetSecondaryElectron = 0.;
  fEnergyExitTargetSecondaryMuonPlus = 0.;
  fEnergyExitTargetSecondaryMuonMinus = 0.;

  fAnglePrimary = 0.;
  fAngleSecondaryPositron = 0.;
  fAngleSecondaryElectron = 0.;
  fAngleSecondaryMuonPlus = 0.;
  fAngleSecondaryMuonMinus = 0.;

  fNbPositronsTarget = fNbElectronsTarget = fNbMuonsPlusTarget = fNbMuonsMinusTarget = 0;
   
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }       

  // do not save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nbEvents = aRun->GetNumberOfEvent();
  if (nbEvents == 0) return;
  
  //G4Material* material = fDetector->GetMaterialT();
  G4double length  = fDetector->GetSizeZW();
   
  G4ParticleDefinition* particle = fPrimary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double ePrimary = fPrimary->GetParticleGun()->GetParticleEnergy();
  
  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";
  G4cout << "\n The run was " << nbEvents << " " << partName << " of "
         << G4BestUnit(ePrimary,"Energy") << " through " 
         << G4BestUnit(length,"Length") << ".";
  G4cout << G4endl;

  G4cout 
    << "\n Number of positrons in target = " 
    << fNbPositronsTarget
    << "\n Number of electrons in target = " 
    << fNbElectronsTarget
    << "\n Number of muons plus in target = " 
    << fNbMuonsPlusTarget
    << "\n Number of muons minus in target = " 
    << fNbMuonsMinusTarget
    << G4endl;
  
  //save histograms 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     

  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }   
    
  if (particle->GetPDGCharge() == 0.) return;
   
  G4cout.precision(5);
  
  G4EmCalculator emCal; 

  G4cout.precision(prec);

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::GetEnergyFromRestrictedRange(G4double range,
            G4ParticleDefinition* particle, G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;
    
  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err  = 1., errmax = 0.00001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;
    Energy -= dE;
    r = emCal.GetRangeFromRestricteDEDX(Energy,particle,material);
    dr = r - range;          
    dEdx = emCal.GetDEDX(Energy,particle,material);
    dE = dEdx*dr;
    err = std::abs(dE)/Energy;    
  }
  if (iter == itermax) {
    G4cout 
    << "\n  ---> warning: RunAction::GetEnergyFromRestRange() did not converge"
    << "   Etry = " << G4BestUnit(Etry,"Energy")
    << "   Energy = " << G4BestUnit(Energy,"Energy")
    << "   err = " << err
    << "   iter = " << iter << G4endl;
  }         
         
  return Energy;         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::GetEnergyFromCSDARange(G4double range,
            G4ParticleDefinition* particle, G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;
    
  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err  = 1., errmax = 0.00001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;
    Energy -= dE;
    r = emCal.GetCSDARange(Energy,particle,material);
    dr = r - range;          
    dEdx = emCal.ComputeTotalDEDX(Energy,particle,material);
    dE = dEdx*dr;
    err = std::abs(dE)/Energy;
  }
  if (iter == itermax) {
    G4cout 
    << "\n  ---> warning: RunAction::GetEnergyFromCSDARange() did not converge"
    << "   Etry = " << G4BestUnit(Etry,"Energy")
    << "   Energy = " << G4BestUnit(Energy,"Energy")
    << "   err = " << err
    << "   iter = " << iter << G4endl;
  }         
         
  return Energy;         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
