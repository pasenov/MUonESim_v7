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
// $Id: DetectorConstruction.hh 66241 2012-12-13 18:34:42Z gunter $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4Cache.hh"
#include "G4ElectricField.hh"
#include "G4EqMagElectricField.hh"


class MagneticField;

class G4LogicalVolume;
class G4Material;
class G4UniformMagField;
class G4UniformElectricField;
class DetectorMessenger;
class ElectricFieldSetup;
class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver;
class FieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
     virtual G4VPhysicalVolume* Construct();
     virtual void ConstructSDandField();
     
     void SetSizeXW     (G4double);
     void SetSizeYW     (G4double);
     void SetSizeZW     (G4double);
     void SetSizeXT     (G4double);
     void SetSizeYT     (G4double);
     void SetSizeZT     (G4double);
     void SetMaterialW (const G4String&);
     void SetMaterialT (const G4String&);
  

     void GetFieldValue (G4double Point[4], G4double* Efield);


     void UpdateGeometry();
     
  public:
  
     const
     G4VPhysicalVolume* GetWorld()      {return fPBoxW;};
     G4VPhysicalVolume* GetTarget()       {return fPBoxT;};
 
     G4LogicalVolume* GetLWorld()      {return fLBoxW;};
     G4LogicalVolume* GetLTarget()       {return fLBoxT;};
                    
                    
     G4double           GetSizeXW()       {return fWorldSizeX;};
     G4double           GetSizeYW()       {return fWorldSizeY;};
     G4double           GetSizeZW()       {return fWorldSizeZ;};

     G4double           GetSizeXT()       {return fTargetSizeX;};
     G4double           GetSizeYT()       {return fTargetSizeY;};
     G4double           GetSizeZT()       {return fTargetSizeZ;};

     G4double           GetElField()   {return fElFieldz;};

     G4Material*        GetMaterialW()   {return fMaterialW;};
     G4Material*        GetMaterialT()   {return fMaterialT;};
     
     void               PrintParameters();
                       
  private:
  
     G4VPhysicalVolume* fPBoxW;
     G4VPhysicalVolume* fPBoxT;

     G4LogicalVolume*   fLBoxW;
     G4LogicalVolume*   fLBoxT;
     
     G4double           zOff; //Offset in z-direction
     G4double           fWorldSizeX;
     G4double           fWorldSizeY;
     G4double           fWorldSizeZ;
     G4double           fTargetSizeX;
     G4double           fTargetSizeY;
     G4double           fTargetSizeZ;

     G4double           fElFieldz;

     G4double		pos_xT, pos_yT, pos_zT;

     G4Material*        fMaterialW;
     G4Material*        fMaterialT;      

     
     DetectorMessenger* fDetectorMessenger;
     G4Cache<ElectricFieldSetup*> fEmFieldSetup;


     G4FieldManager*         fLocalFieldManager;
     G4EqMagElectricField*   fLocalEquation;
     G4ChordFinder*          fLocalChordFinder;
     G4ElectricField*        fElField;
     G4MagIntegratorStepper* fLocalStepper;
     G4MagInt_Driver*        fIntgrDriver;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();
     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

