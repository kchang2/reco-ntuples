#ifndef aobdata_h
#define aobdata_h

#include "TObject.h"


class ARecHit : public TObject
{
public:
  ARecHit() :  layer(0),wafer(0),cell(0),
                x(0.),y(0.),z(0.),
                eta(-1000.),phi(-1000.),pt(-1000.),
                energy(-1000.),time(-1),thickness(0.),
                isHalf(false), flags(0), cluster2d(-1)
  {
  }
 ARecHit(int i_layer, int i_wafer, int i_cell,
	 float i_x, float i_y, float i_z,
     float i_eta, float i_phi, float i_pt,
	 float i_energy, float i_time, float i_thickness,
	 bool i_isHalf, int i_flags, int i_cluster2d):
    layer(i_layer),wafer(i_wafer),cell(i_cell),
    x(i_x),y(i_y),z(i_z),
    eta(i_eta),phi(i_phi),pt(i_pt),
    energy(i_energy),time(i_time),thickness(i_thickness),
    isHalf(i_isHalf),
    flags(i_flags),
    cluster2d(i_cluster2d)
    {
  }


  virtual ~ARecHit() { }

  int   layer,wafer,cell;
  float x,y,z;
  float eta,phi,pt;
  float energy,time,thickness;
  bool  isHalf;
  int   flags;
  int   cluster2d;

  ClassDef(ARecHit,1)
};


typedef std::vector<ARecHit> ARecHitCollection;

#endif
