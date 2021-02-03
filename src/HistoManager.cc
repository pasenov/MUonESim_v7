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
//
// $Id: HistoManager.cc 72242 2013-07-12 08:44:19Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>
#include <sstream>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  :fFileName("MUonESim")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);   //enable inactivation of histograms

  
  // Define histograms start values
  // const G4int kNbofStrips = 1016;
  const G4int kMaxHisto = 33;
  const G4int kMaxHisto2 = 11;

  const G4String id[] = { "000", "001", "002", "003", "004", "005", "006", "007", "008", "009", "010", "011", "012", "013", "014", "015", "016", "017", "018", "019", "020", "021", "022", "023", "024", "025", "026", "027", "028", "029", "030", "031", "032" };
  const G4String id2[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10" };

  const G4String title[] =
                { "dummy",                                                       //0
                  "Energy of primary muons at the exit of the target (simulation)",   //1
                  "Energy of secondary positrons at the exit of the target (simulation)",   //2
                  "Energy of secondary electrons at the exit of the target (simulation)",   //3
                  "Energy of secondary #mu^{+} at the exit of the target (simulation)",   //4
                  "Energy of secondary #mu^{-} at the exit of the target (simulation)",   //5
                  "Angle of primary muons at the exit of the target (simulation, smaller than 30 mrad)",   //6
                  "Angle of secondary positrons at the exit of the target (simulation, smaller than 30 mrad)",   //7
                  "Angle of secondary electrons at the exit of the target (simulation, smaller than 30 mrad)",   //8
                  "Angle of secondary #mu^{+} at the exit of the target (simulation, smaller than 30 mrad)",   //9
                  "Angle of secondary #mu^{-} at the exit of the target (simulation, smaller than 30 mrad)",   //10
                  "Process ID for the creation of secondary particles (simulation)",  //11
                  "Energy of all charged particles at the exit of the target (simulation)",   //12
                  "Angle of all charged particles at the exit of the target (simulation)",   //13
                  "Energy of charged particles at the exit of the target, E > 0.5 GeV (simulation)",   //14
                  "Angle of charged particles at the exit of the target, E > 0.5 GeV (simulation)",   //15
                  "Energy of charged particles at the exit of the target, E > 1 GeV (simulation)",   //16
                  "Angle of charged particles at the exit of the target, E > 1 GeV (simulation)",   //17
                  "Multiplicity of charged tracks with angle < 30 mrad (simulation)",   //18
                  "Multiplicity of charged tracks with angle < 30 mrad, E > 0.5 GeV (simulation)",   //19
                  "Multiplicity of charged tracks with angle < 30 mrad, E > 1 GeV (simulation)",   //20
                  "Energy of secondary e at the exit of the target, ionization (simulation)",   //21   
                  "Energy of secondary e at the exit of the target, pair production (simulation)",   //22
                  "Energy of secondary e at the exit of the target, ionization, E > 0.5 GeV (simulation)",   //23  
                  "Energy of secondary e at the exit of the target, pair production, E > 0.5 GeV (simulation)",   //24 
                  "Energy of secondary e at the exit of the target, ionization, E > 1 GeV (simulation)",   //25 
                  "Energy of secondary e at the exit of the target, pair production, E > 1 GeV (simulation)",   //26
                  "Angle of secondary e at the exit of the target, ionization (simulation)",   //27  
                  "Angle of secondary e at the exit of the target, pair production (simulation)",   //28
                  "Angle of secondary e at the exit of the target, ionization, E > 0.5 GeV (simulation)",   //29 
                  "Angle of secondary e at the exit of the target, pair production, E > 0.5 GeV (simulation)",   //30 
                  "Angle of secondary e at the exit of the target, ionization, E > 1 GeV (simulation)",   //31 
                  "Angle of secondary e at the exit of the target, pair production, E > 1 GeV (simulation)"    //32
                };

  const G4String title2[] =
                { "dummy",                                                     //0
                  "Angle of primary muon vs. angle of secondary positron (simulation)",   //1
                  "Angle of primary muon vs. angle of secondary electron (simulation)",   //2
                  "Angle of primary muon vs. angle of secondary #mu^{+} (simulation)",    //3
                  "Angle of primary muon vs. angle of secondary #mu^{-} (simulation)",    //4
                  "Angle of primary muon vs. angle of secondary e, ionization (simulation)",   //5
                  "Angle of primary muon vs. angle of secondary e, pair production (simulation)",   //6
                  "Angle of primary muon vs. angle of secondary e, ionization, E > 0.5 GeV (simulation)",   //7
                  "Angle of primary muon vs. angle of secondary e, pair production, E > 0.5 GeV (simulation)",   //8
                  "Angle of primary muon vs. angle of secondary e, ionization, E > 1 GeV (simulation)",   //9
                  "Angle of primary muon vs. angle of secondary e, pair production, E > 1 GeV (simulation)"    //10
                };
            
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = -0.5;
  G4double vmax = 100.5;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }

  for (G4int k2=0; k2<kMaxHisto2; k2++) {
    G4int ih2 = analysisManager->CreateH2(id2[k2], title2[k2], nbins, vmin, vmax, nbins, vmin, vmax);
    analysisManager->SetH2Activation(ih2, false);
  }

  G4int ih0 = analysisManager->CreateH2("0", "dummy", nbins, vmin, vmax,  //real point
						      nbins, vmin, vmax);   //reconstructed point
  analysisManager->SetH2Activation(ih0, false);
  
  // nTuples
  analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetFirstNtupleId(1);       
  analysisManager->CreateNtuple("1", "Tuple 1");
  
  std::ostringstream os0;
  os0 <<"Dummy";
  analysisManager->CreateNtupleDColumn(os0.str());
  
  std::ostringstream os1;
  os1 <<"EnergyPrimaryMuonExitTarget";
  analysisManager->CreateNtupleDColumn(os1.str());
  
  std::ostringstream os2;
  os2 <<"EnergySecondaryPositronExitTarget";
  analysisManager->CreateNtupleDColumn(os2.str());
  
  std::ostringstream os3;
  os3 <<"EnergySecondaryElectronExitTarget";
  analysisManager->CreateNtupleDColumn(os3.str());
  
  std::ostringstream os4;
  os4 <<"EnergySecondaryMuonPlusExitTarget";
  analysisManager->CreateNtupleDColumn(os4.str());
  
  std::ostringstream os5;
  os5 <<"EnergySecondaryMuonMinusExitTarget";
  analysisManager->CreateNtupleDColumn(os5.str());
  
  std::ostringstream os6;
  os6 <<"AnglePrimaryMuonExitTarget";
  analysisManager->CreateNtupleDColumn(os6.str());
  
  std::ostringstream os7;
  os7 <<"AngleSecondaryPositronExitTarget";
  analysisManager->CreateNtupleDColumn(os7.str());
  
  std::ostringstream os8;
  os8 <<"AngleSecondaryElectronExitTarget";
  analysisManager->CreateNtupleDColumn(os8.str());
  
  std::ostringstream os9;
  os9 <<"AngleSecondaryMuonPlusExitTarget";
  analysisManager->CreateNtupleDColumn(os9.str());
  
  std::ostringstream os10;
  os10 <<"AngleSecondaryMuonMinusExitTarget";
  analysisManager->CreateNtupleDColumn(os10.str());
  
  std::ostringstream os11;
  os11 <<"ProcessID";
  analysisManager->CreateNtupleDColumn(os11.str());
  
  std::ostringstream os12;
  os12 <<"EnergyChargedExitTarget";
  analysisManager->CreateNtupleDColumn(os12.str());
  
  std::ostringstream os13;
  os13 <<"AngleChargedExitTarget";
  analysisManager->CreateNtupleDColumn(os13.str());
  
  std::ostringstream os14;
  os14 <<"EnergyChargedCut05GeVExitTarget";
  analysisManager->CreateNtupleDColumn(os14.str());
  
  std::ostringstream os15;
  os15 <<"AngleChargedCut05GeVExitTarget";
  analysisManager->CreateNtupleDColumn(os15.str());
  
  std::ostringstream os16;
  os16 <<"EnergyChargedCut1GeVExitTarget";
  analysisManager->CreateNtupleDColumn(os16.str());
  
  std::ostringstream os17;
  os17 <<"AngleChargedCut1GeVExitTarget";
  analysisManager->CreateNtupleDColumn(os17.str());
  
  std::ostringstream os18;
  os18 <<"NbTracks";
  analysisManager->CreateNtupleDColumn(os18.str());
  
  std::ostringstream os19;
  os19 <<"NbTracksCut05GeV";
  analysisManager->CreateNtupleDColumn(os19.str());
  
  std::ostringstream os20;
  os20 <<"NbTracksCut1GeV";
  analysisManager->CreateNtupleDColumn(os20.str());
  
  std::ostringstream os21;
  os21 <<"EnergyIoni";
  analysisManager->CreateNtupleDColumn(os21.str());
  
  std::ostringstream os22;
  os22 <<"EnergyPP";
  analysisManager->CreateNtupleDColumn(os22.str());
  
  std::ostringstream os23;
  os23 <<"EnergyIoniCut05GeV";
  analysisManager->CreateNtupleDColumn(os23.str());
  
  std::ostringstream os24;
  os24 <<"EnergyPPCut05GeV";
  analysisManager->CreateNtupleDColumn(os24.str());
  
  std::ostringstream os25;
  os25 <<"EnergyIoniCut1GeV";
  analysisManager->CreateNtupleDColumn(os25.str());
  
  std::ostringstream os26;
  os26 <<"EnergyPPCut1GeV";
  analysisManager->CreateNtupleDColumn(os26.str());
  
  std::ostringstream os27;
  os27 <<"AngleIoni";
  analysisManager->CreateNtupleDColumn(os27.str());
  
  std::ostringstream os28;
  os28 <<"AnglePP";
  analysisManager->CreateNtupleDColumn(os28.str());
  
  std::ostringstream os29;
  os29 <<"AngleIoniCut05GeV";
  analysisManager->CreateNtupleDColumn(os29.str());
  
  std::ostringstream os30;
  os30 <<"AnglePPCut05GeV";
  analysisManager->CreateNtupleDColumn(os30.str());
  
  std::ostringstream os31;
  os31 <<"AngleIoniCut1GeV";
  analysisManager->CreateNtupleDColumn(os31.str());
  
  std::ostringstream os32;
  os32 <<"AnglePPCut1GeV";
  analysisManager->CreateNtupleDColumn(os32.str());
  
  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
