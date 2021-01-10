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
  const G4int kMaxHisto = 12;
  const G4int kMaxHisto2 = 5;

  const G4String id[] = { "000", "001", "002", "003", "004", "005", "006", "007", "008", "009", "010", "011" };
  const G4String id2[] = { "0", "1", "2", "3", "4" };

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
                  "Process ID for the creation of secondary particles (simulation)"   //11
                };

  const G4String title2[] =
                { "dummy",                                                     //0
                  "Angle of primary muon vs. angle of secondary positron (simulation)",   //1
                  "Angle of primary muon vs. angle of secondary electron (simulation)",   //2
                  "Angle of primary muon vs. angle of secondary #mu^{+} (simulation)",    //3
                  "Angle of primary muon vs. angle of secondary #mu^{-} (simulation)"     //4
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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
