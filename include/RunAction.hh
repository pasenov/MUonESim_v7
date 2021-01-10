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
/// \file electromagnetic/TestEm18/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class G4ParticleDefinition;
class G4Material;

class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddREnergyExitTargetPrimary (G4double edep)
                   {fEnergyExitTargetPrimary += edep;};

    void AddREnergyExitTargetSecondaryPositron (G4double edep)
                   {fEnergyExitTargetSecondaryPositron += edep;};

    void AddREnergyExitTargetSecondaryElectron (G4double edep)
                   {fEnergyExitTargetSecondaryElectron += edep;};

    void AddREnergyExitTargetSecondaryMuonPlus (G4double edep)
                   {fEnergyExitTargetSecondaryMuonPlus += edep;};

    void AddREnergyExitTargetSecondaryMuonMinus (G4double edep)
                   {fEnergyExitTargetSecondaryMuonMinus += edep;};

    void AddRAnglePrimary (G4double sangle)
                   {fAnglePrimary += sangle;};

    void AddRAngleSecondaryPositron (G4double sangle)
                   {fAngleSecondaryPositron += sangle;};

    void AddRAngleSecondaryElectron (G4double sangle)
                   {fAngleSecondaryElectron += sangle;};

    void AddRAngleSecondaryMuonPlus (G4double sangle)
                   {fAngleSecondaryMuonPlus += sangle;};

    void AddRAngleSecondaryMuonMinus (G4double sangle)
                   {fAngleSecondaryMuonMinus += sangle;};

    void AddNbPositronsTarget ()
                 {fNbPositronsTarget++;};

    void AddNbElectronsTarget ()
                 {fNbElectronsTarget++;};

    void AddNbMuonsPlusTarget ()
                 {fNbMuonsPlusTarget++;};

    void AddNbMuonsMinusTarget ()
                 {fNbMuonsMinusTarget++;};
 
               
  public:
    G4double GetEnergyFromRestrictedRange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);
                       
    G4double GetEnergyFromCSDARange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);                 
                 
  private:
    G4double fEnergyExitTargetPrimary, fEnergyExitTargetSecondaryPositron, fEnergyExitTargetSecondaryElectron, fEnergyExitTargetSecondaryMuonPlus, fEnergyExitTargetSecondaryMuonMinus;
    G4double fAnglePrimary, fAngleSecondaryPositron, fAngleSecondaryElectron, fAngleSecondaryMuonPlus, fAngleSecondaryMuonMinus;

    G4int fNbPositronsTarget, fNbElectronsTarget, fNbMuonsPlusTarget, fNbMuonsMinusTarget;

    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    HistoManager*           fHistoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

