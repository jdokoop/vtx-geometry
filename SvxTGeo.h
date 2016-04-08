// $Id: SvxTGeo.h,v 1.4 2014/05/06 17:26:21 adare Exp $

#ifndef __SVXTGEO_H__
#define __SVXTGEO_H__

#include <vector>
#include <TMatrixD.h>

class TGeoMaterial;
class TGeoMedium;
class TGeoVolume;
class TGeoNode;
class TGeoManager;
class TGeoTranslation;
class TGeoCombiTrans;
class TPolyLine;

using namespace std;

class SvxTGeo
{
public:
  SvxTGeo();
  ~SvxTGeo();
  TGeoVolume *MakeTopVolume(double x=100, double y=100, double z=100 /* cm */);
  TGeoVolume *MakeBox(double xhw, double yhw, double zhw, const char *name = 0);
  void AddVolume(TGeoVolume *parent, TGeoVolume *daughter,
                 double x, double y, double z,
                 double phi, double theta, double psi);
  void AddSensor(int lyr, int ldr, int sns);
  void AddSensors();
  void InitMaterials();
  void ReadParFile(const char *filename);
  void WriteParFile(const char *filename);
  void WritePar(ofstream &fs, const int par);
  void WritePar(ofstream &fs, const float par);
  void WritePars(ofstream &fs, vector<int> &vec);
  void WritePars(ofstream &fs, vector<float> &vec);
  void WriteGlobalPars(ofstream &fs);
  void WriteBarrelPars(ofstream &fs);
  void WriteSensorPars(ofstream &fs);
  void TranslateLadder(int layer, int ladder, float x, float y, float z);
  void MoveLadderRadially(int layer, int ladder, float dr /*cm*/);
  void RotateLadder(int layer, int ladder, float polarx, float polary, float phi);
  void RotateLadderRPhi(int layer, int ladder, float rphi);
  void GetSensorXYZ(int lyr, int ldr, int sns, double *xyz);
  float LayerRadius(int layer);
  float SensorRadius(int layer, int ladder, int sensor);
  float SensorPhiRad(int layer, int ladder, int sensor);
  float SensorPhiDeg(int layer, int ladder, int sensor);
  float GetLadderPhiTilt(int layer, int ladder); // In radians
  TGeoNode *SensorNode(int lyr, int ldr, int sns);
  TPolyLine *LadderOutlineXY(int lyr, int ldr);
  void AddLadder(int lyr, int ldr,
                 double x, double y, double zoffset,
                 double phi, double theta=0, double psi=0);

  TGeoManager *GeoManager()
  {
    return fGeoMgr;
  }

  void
  GetSensorHalfWidths(const int layer, double &xhw, double &yhw, double &zhw)
  {
    xhw = fSensorXHW[layer];
    yhw = fSensorYHW[layer];
    zhw = fSensorZHW[layer];
  }

  int GetNLayers()
  {
    return fNLayers;
  }

  int GetNLadders(int lyr)
  {
    return fNLadders.at(lyr);
  }

  int GetNSensors(int lyr)
  {
    return fNSensors.at(lyr);
  }

  void SetVerbosity(int v)
  {
    fVerbosity = v;
  }

  void SetPisaFileColumnWidth(int w)
  {
    fColWidth = w;
  }

  void SetPisaFilePrecision(int p)
  {
    fPrec = p;
  }

  struct GBox
  {
    GBox() :
      x(0.), y(0.), z(0.),
      phi(0.), theta(0.), psi(0.),
      xhw(0.), yhw(0.), zhw(0.)
    {
      R.ResizeTo(3,3);
    }
    ~GBox() {};

    double x,y,z;           // Position at centroid (relative to parent node)
    double phi, theta, psi; // Euler angles [deg]
    double xhw, yhw, zhw;   // Half-widths [cm]
    TMatrixD R;             // 3x3 rotation matrix Rx * Ry * Rz.
  };

  std::vector<GBox> sensors;
  int indx[8][48][12]; // Allocate space for 2x max to accommodate new objects

protected:
  TGeoManager  *fGeoMgr;
  TGeoVolume   *fTopVolume;
  TGeoMaterial *fVacuumMaterial;
  TGeoMaterial *fSiliconMaterial;
  TGeoMaterial *fAluminumMaterial;
  TGeoMedium   *fVacuumMedia;
  TGeoMedium   *fSiliconMedia;
  TGeoMedium   *fAluminumMedia;
  bool fNewGeo;
  int fComponentId; // Unique identifier for all elements
  int fVerbosity;
  int fColWidth; // Column width of fields in PISA output file. Default 12
  int fPrec;     // Precision of fields in PISA output file. Default 5 (0.1 um)

  // SVX parameters
  std::vector<float> fCage1;       // 06: ? Not used.
  std::vector<float> fCage2;       // 07: ? Not used.
  int fNLayers;                    // 10: # layers.
  std::vector<float> fRadii;       // 11: Layer radii, stagger.
  std::vector<int>   fNSensors;    // 13: # sensors/ladder in layer i.
  std::vector<float> fSensorXHW;   // 14: Sensor x half-width in layer i
  std::vector<float> fSensorYHW;   // 15: Sensor y half-width in layer i
  std::vector<float> fSensorZHW;   // 16: Sensor z half-width in layer i
  std::vector<float> fSensorZGap;  // 17: Gap btwn. sensors in layer i
  std::vector<float> fX0add;       // 18: ? Not used.
  std::vector<float> fDPhi;        // 19: ? Not used.
  std::vector<float> fTilt;        // 20: ? Not used.
  std::vector<int>   fNArms;       // 21, 24, 27, 30: Number of arms in layers 0-3.
  std::vector<int>   fNLaddersE;   // 22, 25, 28, 31: ladders in E,W (or W,E?) arm.
  std::vector<int>   fNLaddersW;
  std::vector<int>   fNLadders;

  ClassDef(SvxTGeo,1)
};

// ----------------------------------------------------------------------------

#endif // __SVXTGEO_H__












