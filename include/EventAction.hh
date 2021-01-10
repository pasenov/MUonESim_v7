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
/// \file electromagnetic/TestEm18/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*, RunAction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void   EndOfEventAction(const G4Event*);

    void AddEnergyExitTargetPrimary(G4double edep)            {fEnergyExitTargetPrimary  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryPositron(G4double edep)  {fEnergyExitTargetSecondaryPositron  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryElectron(G4double edep)  {fEnergyExitTargetSecondaryElectron  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryMuonPlus(G4double edep)  {fEnergyExitTargetSecondaryMuonPlus  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryMuonMinus(G4double edep) {fEnergyExitTargetSecondaryMuonMinus  += edep;}; // If angle < 30 mrad

    void AddAnglePrimary(G4double theta)            	      {if (fBoolPrim == 0)   { fAnglePrimary  += theta; fBoolPrim = 1; } }; // If angle < 30 mrad
    void AddAngleSecondaryPositron(G4double theta)            {if (fBoolSecPositron == 0)   {fAngleSecondaryPositron  += theta; fBoolSecPositron = 1;} }; // If angle < 30 mrad
    void AddAngleSecondaryElectron(G4double theta)            {if (fBoolSecElectron == 0)   {fAngleSecondaryElectron  += theta; fBoolSecElectron = 1;} }; // If angle < 30 mrad
    void AddAngleSecondaryMuonPlus(G4double theta)            {if (fBoolSecMuonPlus == 0)   {fAngleSecondaryMuonPlus  += theta; fBoolSecMuonPlus = 1;} }; // If angle < 30 mrad
    void AddAngleSecondaryMuonMinus(G4double theta)           {if (fBoolSecMuonMinus == 0)  {fAngleSecondaryMuonMinus += theta; fBoolSecMuonMinus = 1;} }; // If angle < 30 mrad

    /*void AddMomentumDirectionPrimExT(G4double dirx, G4double diry, G4double dirz)  {fMomDirPx += dirx; fMomDirPy += diry; fMomDirPz += dirz;};
    void AddMomentumDirectionSecExT(G4double dirx, G4double diry, G4double dirz)  {fMomDirSx += dirx; fMomDirSy += diry; fMomDirSz += dirz;};*/
       
  private:
    DetectorConstruction* fDetectorconstruction;
    RunAction*    fRunAction;
    
    G4double      fEnergyExitTargetPrimary;
    G4double      fEnergyExitTargetSecondaryPositron;
    G4double      fEnergyExitTargetSecondaryElectron;
    G4double      fEnergyExitTargetSecondaryMuonPlus;
    G4double      fEnergyExitTargetSecondaryMuonMinus;

    G4double      fAnglePrimary;
    G4double      fAngleSecondaryPositron;
    G4double      fAngleSecondaryElectron;
    G4double      fAngleSecondaryMuonPlus;
    G4double      fAngleSecondaryMuonMinus;

    G4int         fBoolPrim, fBoolSecPositron, fBoolSecElectron, fBoolSecMuonPlus, fBoolSecMuonMinus; //a boolean which becomes 1 when a primary or secondary, respectively, is found with angle < 30 mrad

    /*G4double      fMomDirPx, fMomDirPy, fMomDirPz;
    G4double      fMomDirSx, fMomDirSy, fMomDirSz;*/
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
