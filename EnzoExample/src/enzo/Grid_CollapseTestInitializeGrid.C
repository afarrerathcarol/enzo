/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COLLAPSE TEST)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define NTHETA 1000
#define NR 1000

/********************* PROTOTYPES *********************/

int GetUnits(float *DensityUnits, float *LengthUnits,
     float *TemperatureUnits, float *TimeUnits,
     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
float gasdev();

// Used to compute Bonner-Ebert density profile
double BE(double r);
double q(double r);
double Ang(double a1, double a2, double R, double r);

// Returns random velocity from Maxwellian distribution
double Maxwellian(double c_tilda, double vel_unit, double mu, double gamma);
double ERF(double x);

int ComputeRadialVelocity(float density, double mass, float r_init, 
      double VelocityUnits, double LengthUnits,
      float dm_velrot_corr, double radius_vr[], double Vr[], 
      double exterior_rho[], int Npts);

double ComputeBubbleVelocityInducedByDarkMatter(const float *PointSourceGravityPosition, 
    const float PointSourceGravityCoreRadius, const float *center_of_bubble, 
    float *BubbleVelocityTowardCenter, const float DensityUnits, const float LengthUnits, 
    const float TemperatureUnits, const float TimeUnits, const float VelocityUnits, const float SolarMass,
    const float GravConst, const float SphereDensity, const int MetallicityDistributionCase, const float SphereCoreRadius);


/*******************************************************/

static int CollapseTestParticleCount = 0;

static float CosmologySimulationOmegaBaryonNow       = 0.0463;
static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;
// static float CollapseTestInitialMetallicity = 1e-6;
static float CollapseTestInitialMetallicity = 1e-20;
static float SphereColour = 1e-6;

int grid::CollapseTestInitializeGrid(
          const int NumberOfSpheres,
          const FLOAT SphereRadius[MAX_SPHERES],
          const FLOAT SphereCoreRadius[MAX_SPHERES],
          const float SphereDensity[MAX_SPHERES],
          const float SphereTemperature[MAX_SPHERES],
          const float SphereMetallicity[MAX_SPHERES],
          const FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
          const float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
           float SphereFracKeplerianRot[MAX_SPHERES],
          const float SphereTurbulence[MAX_SPHERES],
          const float SphereDispersion[MAX_SPHERES],
          const float SphereCutOff[MAX_SPHERES],
          const float SphereAng1[MAX_SPHERES],
          const float SphereAng2[MAX_SPHERES],
          const int   SphereNumShells[MAX_SPHERES],
          const int   SphereType[MAX_SPHERES],
          const int   SphereSubType[MAX_SPHERES],
          const float SphereSubTypeSlope[MAX_SPHERES],
          const int   SphereConstantPressure[MAX_SPHERES],
          const int   SphereSmoothSurface[MAX_SPHERES],
          const float SphereSmoothRadius[MAX_SPHERES],
          const float SphereHII[MAX_SPHERES],
          const float SphereHeII[MAX_SPHERES],
          const float SphereHeIII[MAX_SPHERES],
          const float SphereH2I[MAX_SPHERES],
          const int   SphereUseParticles,
           float ParticleMeanDensity,
          const float UniformVelocity[MAX_DIMENSION],
          const int   SphereUseColour,
          const int   SphereUseMetals,
          const float InitialTemperature, 
          const float InitialDensity, 
          const int   level,
          const float CollapseTestInitialFractionHII, 
          const float CollapseTestInitialFractionHeII,
          const float CollapseTestInitialFractionHeIII, 
          const float CollapseTestInitialFractionHM,
          const float CollapseTestInitialFractionH2I, 
          const float CollapseTestInitialFractionH2II, 
#ifdef TRANSFER
           float RadHydroRadiation,
#endif          
           float center_of_bubble[MAX_SPHERES][MAX_BUBBLES][MAX_DIMENSION], 
          const int   MetallicityDistributionCase[MAX_SPHERES], 
          const int   NumberOfBubbles[MAX_SPHERES],
          const float BubbleVelocity[MAX_SPHERES][MAX_BUBBLES][MAX_DIMENSION],
          const float BubbleDensity[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleMetallicity[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleTemperature[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleRadius[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleInnerRadius[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleInnerDensity[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleOuterDensity[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleInnerMetallicity[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleOuterMetallicity[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleInnerTemperature[MAX_SPHERES][MAX_BUBBLES],
          const float BubbleOuterTemperature[MAX_SPHERES][MAX_BUBBLES], 
          const int   VelocityTowardCenter[MAX_SPHERES])
{
  /* declarations */
  int dim, i, j, k, m, field, sphere, size;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum, DINum, DIINum, HDINum, MetalNum, MetalIaNum;
  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, PhotoGammaNum;
#ifdef TRANSFER
  int EgNum;
#endif  
  float xdist,ydist,zdist;

  /* create fields */
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int ivel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
#ifdef TRANSFER
  if (RadiativeTransferFLD > 1) {
    FieldType[EgNum = NumberOfBaryonFields++] = RadiationFreq0;
    if (RadiativeCooling) {
      FieldType[ kphHINum = NumberOfBaryonFields++] = kphHI;
      FieldType[ PhotoGammaNum = NumberOfBaryonFields++] = PhotoGamma;
      if (RadiativeTransferHydrogenOnly == FALSE) {
        FieldType[ kphHeINum = NumberOfBaryonFields++] = kphHeI;
        FieldType[ kphHeIINum = NumberOfBaryonFields++] = kphHeII;
      }
      if (MultiSpecies > 1)
        FieldType[ kdissH2INum = NumberOfBaryonFields++] = kdissH2I;
     }
  }
#endif


  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }
  if (SphereUseMetals) {
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity; /* 2018.02.12 It was originally SNColour and I changed it to Metallicity. */
  /* 2018.03.06 Changed it back to SNColour */
  // printf("Metallicity: %d\n", Metallicity);
  }
  if (StarMakerTypeIaSNe) {
    FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity; /* 2018.03.06 Also modified the code in Grid_SolveRateAndCoolEquations.C
   to be consistent.*/
  }

  int ColourNum = NumberOfBaryonFields;
  if (SphereUseColour) {
    FieldType[NumberOfBaryonFields++] = SNColour; /* fake it with metals, 2018.03.05 be careful when use color field. */
    printf("UseColour\n");
  }

  /* Return if this doesn't concern us. */
  if (ProcessorNumber != MyProcessorNumber) {
    NumberOfParticles = (SphereUseParticles > 0) ? 1 : 0; // if SphereUseParticles > 0 => NumberOfParticles = 1, else NumberOfParticles = 0
  for (dim = 0; dim < GridRank; dim++)
    NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
  }

  /* Set various units. */

  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8,
     pi = 3.14159, mh = 1.67e-24, kboltz = 1.381e-16, LightSpeed = 2.9979e10;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, mu = 0.6;

  FLOAT a, dadt, ExpansionFactor = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CriticalDensity = 2.78e11*pow(HubbleConstantNow, 2); // in Msolar/Mpc^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
  } else {
    CriticalDensity = 2.78e11*pow(0.74,2); // in Msolar/Mpc^3 for h=0.74
    BoxLength = LengthUnits / 3.086e24;
    HubbleConstantNow = 1.0;
    OmegaMatterNow = 1.0;
  }

  /* Compute NFW profile-info. The parameters are SphereCoreRadius (which is
   the "knee radius") and SphereDensity (the overdensity ar the knee
   radius).  The pressure is computed by integrated the hydro-static
   equation and from this comes the temperature and the dm velocity
   dispersion. */

  #define NFW_POINTS 500
  float NFWMass[NFW_POINTS], NFWPressure[NFW_POINTS], NFWTemp[NFW_POINTS], x1, NFWDensity[NFW_POINTS], NFWSigma[NFW_POINTS], m200;
  FLOAT NFWRadius[NFW_POINTS];
  double dpdr = 0, dpdr_old;
  sphere = 0;
  m200   = 0;
  NFWPressure[0] = 1.0 * kboltz * InitialTemperature / (mu * mh);
  FILE *fptr = fopen("NFWProfile.out", "w");
  for (i = 0; i < NFW_POINTS; i++) {
  NFWRadius[i] = SphereRadius[sphere]*pow(10, -3*(float(i)/NFW_POINTS));
  x1 = NFWRadius[i]/SphereCoreRadius[sphere];
  NFWDensity[i] = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));
  NFWMass[i] = 4.0*pi*SphereDensity[sphere] * (CriticalDensity/pow(ExpansionFactor, 3)) 
       * pow(SphereCoreRadius[sphere]*BoxLength, 3) * (log(1.0+x1) - x1/(x1+1.0));  // in Msolar
  dpdr_old = dpdr;
  dpdr = GravConst * NFWMass[i] * SolarMass * NFWDensity[i] / pow(NFWRadius[i]*BoxLength*Mpc, 2);
  if (i > 0)
    NFWPressure[i] = NFWPressure[i-1] - 0.5*(dpdr+dpdr_old)*(NFWRadius[i]-NFWRadius[i-1])*BoxLength*Mpc;
  NFWTemp[i]  = NFWPressure[i]*mu*mh/(kboltz*NFWDensity[i]); // in K
  NFWSigma[i] = sqrt(kboltz * NFWTemp[i] / (mu * mh));  // in cm/s
  float mean_overdensity = 3.0*SphereDensity[sphere] / (x1*x1*x1) * (log(1.0+x1) - x1/(x1+1.0));
  fprintf(fptr, "%"ISYM" %"GOUTSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", i, NFWRadius[i], 
    NFWDensity[i], NFWMass[i], NFWPressure[i], NFWTemp[i], NFWSigma[i], mean_overdensity);
  if (mean_overdensity > 200 && m200 == 0)
    m200 = NFWMass[i];
  }
  fprintf(fptr, "#m200 = %"GSYM"\n", m200);
  fclose(fptr);

  /* Set up the velocity of bubbles toward the center (mimic the effect of dark matter halo) 
   Use the cooridinates of the center of bubble as reference, 
   ignore the fact that bubbles have a physical size.
  */

    float BubbleVelocityTowardCenter[MAX_SPHERES][MAX_BUBBLES][MAX_DIMENSION];
    double v_mag = 0.0;
    for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
      for (int BubbleCount = 0; BubbleCount < NumberOfBubbles[sphere]; BubbleCount++) {
        for (dim = 0; dim < MAX_DIMENSION; dim++) {
          BubbleVelocityTowardCenter[sphere][BubbleCount][dim] = 0.0;        
        }
      }
    }

  
    float dist_to_halocenter, xx, yy, zz, norm_factor;
    for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
      if (VelocityTowardCenter[sphere] > 0) {
        for (int BubbleCount = 0; BubbleCount < NumberOfBubbles[sphere]; BubbleCount++) {
          // printf("PointSourceGravityPosition[0]:%f \n", PointSourceGravityPosition[0]);
          // printf("center_of_bubble[0]: %f\n", center_of_bubble[0][0][0]);
          xx = PointSourceGravityPosition[0] - center_of_bubble[sphere][BubbleCount][0]; // code unit
          yy = PointSourceGravityPosition[1] - center_of_bubble[sphere][BubbleCount][1]; // code unit
          zz = PointSourceGravityPosition[2] - center_of_bubble[sphere][BubbleCount][2]; // code unit
          xx *= LengthUnits;
          yy *= LengthUnits;
          zz *= LengthUnits;
          norm_factor = sqrt(xx*xx + yy*yy + zz*zz);
          
          //  printf("norm_factor: %.3e\n", norm_factor);
          //  printf("xx/norm_factor: %.3e\n", xx/norm_factor);
          //  printf("yy/norm_factor: %.3e\n", yy/norm_factor);
          //  printf("zz/norm_factor: %.3e\n", zz/norm_factor);
          //  printf("test:%.3e\n", xx/norm_factor*xx/norm_factor + yy/norm_factor*yy/norm_factor + zz/norm_factor*zz/norm_factor);

          // call function to calculate the velocity toward the center
          v_mag = ComputeBubbleVelocityInducedByDarkMatter(PointSourceGravityPosition, PointSourceGravityCoreRadius, 
            center_of_bubble[sphere][BubbleCount], BubbleVelocityTowardCenter[sphere][BubbleCount],
            DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, SolarMass, 
            GravConst, SphereDensity[sphere], MetallicityDistributionCase[sphere], SphereCoreRadius[sphere]);      
            // printf("v_mag: %.3e\n", v_mag);
          BubbleVelocityTowardCenter[sphere][BubbleCount][0] = VelocityTowardCenter[sphere] * v_mag * xx / norm_factor;
          BubbleVelocityTowardCenter[sphere][BubbleCount][1] = VelocityTowardCenter[sphere] * v_mag * yy / norm_factor;
          BubbleVelocityTowardCenter[sphere][BubbleCount][2] = VelocityTowardCenter[sphere] * v_mag * zz / norm_factor;

          // printf("vx: %.3e, vy: %.3e, vz: %.3e\n", BubbleVelocityTowardCenter[sphere][BubbleCount][0], 
                  // BubbleVelocityTowardCenter[sphere][BubbleCount][1], BubbleVelocityTowardCenter[sphere][BubbleCount][2]);
        }
      }
    }

    

  /* Loop over the set-up twice, once to count the particles, the second
   time to initialize them. */
  int SetupLoopCount, npart = 0;
  for (SetupLoopCount = 0; SetupLoopCount < 1+min(SphereUseParticles, 1); SetupLoopCount++) {

  /* Set densities */
  float BaryonMeanDensity, ParticleCount = 0;
  switch (SphereUseParticles) {
    case 1:
    BaryonMeanDensity = CosmologySimulationOmegaBaryonNow / OmegaMatterNow;
    break;
    case 2:
    BaryonMeanDensity = 1 - CosmologySimulationOmegaBaryonNow / OmegaMatterNow;
    break;
    default:
    BaryonMeanDensity = 1.0; // 2018.02.02 for the background I suppose?
  } // ENDSWITCH SphereUseParticles

  if (!ComovingCoordinates) BaryonMeanDensity = 1.0;

  if (ParticleMeanDensity == FLOAT_UNDEFINED) ParticleMeanDensity = 1.0 - BaryonMeanDensity;
  else BaryonMeanDensity = 1.0 - ParticleMeanDensity;

  /* Set particles. */
  if (SphereUseParticles > 0 && SetupLoopCount > 0) {

    /* If particles already exist (coarse particles), then delete. */
    if (NumberOfParticles > 0) this->DeleteParticles();

    /* Use count from previous loop to set particle number. */
    NumberOfParticles = npart;
    npart = 0;

    /* Allocate space. */
    this->AllocateNewParticles(NumberOfParticles);

    /* Particle values will be set below. */
  } // end: particle initialization

  /* Set up the baryon field. */

  /* compute size of fields */
  size = 1;
  for (dim = 0; dim < GridRank; dim++) size *= GridDimension[dim];

  /* allocate fields */
  if (SetupLoopCount == 0)
    for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];

  /* Loop over the mesh. */
  float density, dens1=0., old_density, Velocity[MAX_DIMENSION], 
    DiskVelocity[MAX_DIMENSION], temperature, temp1, sigma, sigma1, 
    colour, weight, a, DMVelocity[MAX_DIMENSION], metallicity, outer_radius,
    metallicity_to_initialize, temperature_to_intialize, density_to_initialize;

  // /* Added 2018.03.05 */
  // float _to_initialize[2];
  // _to_initialize[0] = 1.0e-20; // metallicity_to_initialize
  // _to_initialize[1] = 10000000; // temperature_to_intialize
  // /* Added 2018.03.05 */

  FLOAT r, rcyl, x, y = 0, z = 0;
  int n = 0, ibin;

  double SphereRotationalPeriod[MAX_SPHERES];
  double VelocityKep = 0;
  double RotVelocity[MAX_DIMENSION];
  float  DMRotVelocityCorrection = 1, GasRotVelocityCorrection = 1;
  double VelocitySound[MAX_SPHERES];
  double SphereMass, SphereCoreMass, SphereCoreDens;
  double alpha, beta, theta;
  double SphereCritMass;
  double Scale_Factor[MAX_SPHERES];
  double term1, term2;
  double radius_vr[NR], vr[NR], exterior_rho[NR], radial_velocity;
  float sin_deltaDisk[MAX_SPHERES];
  double SchwarzschildRadius, CavityRadius, InnerDensity, InnerTemperature,
    ThickenTransitionRadius, BHMass, ScaleHeight, InnerScaleHeight;
  double MidplaneDensity, MidplaneTemperature;
  float HII_Fraction, HeII_Fraction, HeIII_Fraction, H2I_Fraction;



#ifdef TRANSFER
  // if using FLD-based radiation energy density, set the field
  if ((RadiativeTransferFLD > 1)) {
    float RadScaled = RadHydroRadiation/DensityUnits/VelocityUnits/VelocityUnits;
    for (i=0; i<size; i++)
      BaryonField[EgNum][i] = RadScaled;
  }
#endif



  /* Pre-compute cloud properties before looping over mesh */
  for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
    // printf("SphereRadius:%f\n", SphereRadius[sphere]);
    Scale_Factor[sphere] = SphereCutOff[sphere] / SphereRadius[sphere];
    sin_deltaDisk[sphere] = sin(pi * SphereCutOff[sphere] / 180.0);

    if (SphereFracKeplerianRot[sphere] <= 0.0 && 
    (SphereType[sphere] == 7 || SphereType[sphere] == 9) && 
    PointSourceGravityConstant > 0)
    
    SphereFracKeplerianRot[sphere] = 1;

    switch (SphereType[sphere]) {

    case 1:
      SphereMass = (4*pi/3)*pow((SphereRadius[sphere]*LengthUnits), 3) * (SphereDensity[sphere]*DensityUnits);
      printf("mass = %"GSYM", lunit = %"GSYM", dunit = %"GSYM", rho = %"GSYM", r = %"GSYM"\n",
        SphereMass, LengthUnits, DensityUnits, SphereDensity[sphere], SphereRadius[sphere]);
      break;

    case 3:
      SphereMass = m200*SolarMass;
      break;

    case 4:
      // Integral of a Gaussian
      term1 = 0.25 * sqrt(pi) * 0.3536 *  // 2^(-3/2) = 0.3536
        pow((SphereCoreRadius[sphere]), 3) *
        ERF(SphereRadius[sphere]/SphereCoreRadius[sphere]);
      term2 = 0.25 * SphereRadius[sphere] *
        pow(SphereCoreRadius[sphere], 2) *
        PEXP(-0.5 * pow((SphereRadius[sphere]/SphereCoreRadius[sphere]), 2));
      SphereMass = (4*pi*SphereDensity[sphere]*DensityUnits) * pow(LengthUnits, 3) * (term1 - term2);
      break;

    case 5:
      SphereCoreDens = (SphereDensity[sphere]*DensityUnits) * 
          pow(SphereCoreRadius[sphere] / SphereRadius[sphere], -2);
      SphereCoreMass = (4*pi/3) * 
          pow((SphereCoreRadius[sphere]*LengthUnits),3) * (SphereCoreDens);
      SphereMass     = (4*pi*SphereDensity[sphere]*DensityUnits * pow(SphereRadius[sphere]*LengthUnits, 2) *
          ((SphereRadius[sphere] - SphereCoreRadius[sphere]) * LengthUnits)) + SphereCoreMass;
      break;

    case 6:
      SphereMass = pow(VelocitySound[sphere],3) / (sqrt(4*pi*pow((GravConst)/(4*pi),3) * 
        (SphereDensity[sphere]))) * q(SphereCutOff[sphere]);
      SphereMass = SphereMass*DensityUnits*pow(LengthUnits,3);
      break;

    case 7:
      printf("PointSourceGravityConstant[P%"ISYM"][0] = %"GSYM"\n", MyProcessorNumber, PointSourceGravityConstant);
      if (level == 0)
      SphereMass = PointSourceGravityConstant * SolarMass;
      else
      SphereMass = PointSourceGravityConstant * (LengthUnits*VelocityUnits*VelocityUnits) / GravConst;
      // Convert to code units (accel = vel^2 = 1 at dr = 1)
      PointSourceGravityConstant = GravConst * SphereMass / (LengthUnits*VelocityUnits*VelocityUnits);
      printf("PointSourceGravityConstant[P%"ISYM"][1] = %"GSYM"\n", MyProcessorNumber, PointSourceGravityConstant);
      break;

    case 8:
      SphereMass = (4*pi/3)*pow((SphereRadius[sphere]*LengthUnits), 3) * (SphereDensity[sphere]*DensityUnits);
      ComputeRadialVelocity(SphereDensity[sphere], SphereMass, SphereRadius[sphere], VelocityUnits,
            LengthUnits, DMRotVelocityCorrection, radius_vr, vr, exterior_rho, NR);
      GasRotVelocityCorrection = DMRotVelocityCorrection;
      break;

    // Circumbinary BH accretion disk (Lippai et al. 2008)
    case 9:
      printf("PointSourceGravityConstant[P%"ISYM"][0] = %"GSYM"\n", MyProcessorNumber, PointSourceGravityConstant);
      if (level == 0)
      SphereMass = PointSourceGravityConstant * SolarMass;
      else
      SphereMass = PointSourceGravityConstant * (LengthUnits*VelocityUnits*VelocityUnits) / GravConst;
      // Convert to code units (accel = vel^2 = 1 at dr = 1)
      PointSourceGravityConstant = GravConst * SphereMass / (LengthUnits*VelocityUnits*VelocityUnits);
      printf("PointSourceGravityConstant[P%"ISYM"][1] = %"GSYM"\n", MyProcessorNumber, PointSourceGravityConstant);

      BHMass = SphereMass / SolarMass;  // in solar masses
      SchwarzschildRadius = 2.0 * GravConst * SphereMass / (LightSpeed*LightSpeed);
      CavityRadius = 117.0 * SchwarzschildRadius * pow((BHMass/1e6), 0.08);
      InnerDensity = 4.31e-10 * pow((BHMass/1e8), -0.8) *
          pow((CavityRadius/SchwarzschildRadius/1e3), -0.6) / DensityUnits;
      InnerTemperature = 1.7e6 * pow((BHMass/1e6), -0.28);
      InnerScaleHeight = 0.46 * CavityRadius * pow(BHMass/1e6, -0.12) / LengthUnits;
      ThickenTransitionRadius = 1.9e3 * SchwarzschildRadius * pow((BHMass/1e6), 2.0/21) / LengthUnits;
      CavityRadius /= LengthUnits;
      printf("cgs: %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", SchwarzschildRadius, 
        CavityRadius*LengthUnits, InnerDensity*DensityUnits, 
        InnerTemperature, ThickenTransitionRadius*LengthUnits,
        InnerScaleHeight*LengthUnits);
      printf("code: %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", SchwarzschildRadius/LengthUnits, 
        CavityRadius, InnerDensity, InnerTemperature, 
        ThickenTransitionRadius, InnerScaleHeight);
      break;
    

    case 11:
      SphereMass = (4*pi/3)*pow((SphereRadius[sphere]*LengthUnits), 3) * (SphereDensity[sphere]*DensityUnits);
      float BubbleTotalMass = 0.0;
      // printf("SphereMass: %"GSYM"\n", SphereMass);
      for (int bubble=0; bubble < NumberOfBubbles[sphere]; bubble++) {
        if (MetallicityDistributionCase[sphere] == 1) {
          BubbleTotalMass += (4*pi/3)*pow((0.1*SphereRadius[sphere]*LengthUnits), 3) * (BubbleDensity[sphere][bubble]*DensityUnits);              
        } else if (MetallicityDistributionCase[sphere] == 2) {
          BubbleTotalMass += (4*pi/3)*pow((0.1*SphereRadius[sphere]*LengthUnits), 3) * (BubbleDensity[sphere][bubble]*DensityUnits);             
        } else if (MetallicityDistributionCase[sphere] == 3) {
          BubbleTotalMass += (4*pi/3)*pow((0.05*SphereRadius[sphere]*LengthUnits), 3) * (BubbleDensity[sphere][bubble]*DensityUnits);
        } else if (MetallicityDistributionCase[sphere] >= 4) {

          float BaryonMassOfBigSphere, rho_, SphereCoreBaryonMass;
          SphereCoreBaryonMass = 4. * 3.14159 * pow(SphereCoreRadius[sphere]*LengthUnits, 3) / 3. * SphereDensity[sphere] * DensityUnits;
          BaryonMassOfBigSphere = SphereCoreBaryonMass;

          float r_prime = SphereCoreRadius[sphere]*LengthUnits;      
          float deltar = (SphereRadius[sphere]*LengthUnits - SphereCoreRadius[sphere]*LengthUnits) / 1000000.;

          for (int iteration = 0; iteration < 1000000; iteration ++ ) {
            r_prime += deltar;
            BaryonMassOfBigSphere += 4 * 3.14159 * (r_prime * r_prime) * deltar * DensityUnits
                                       * SphereDensity[sphere] * pow((r_prime/(SphereCoreRadius[sphere]*LengthUnits)), -2.2);
            // if (iteration % 100000 == 0) {
              // printf("r_prime: %.3e pc\n", r_prime/3.0857e+18);
              // printf("BaryonMass: %.6e M_sun\n", BaryonMassOfBigSphere/SolarMass);
              // printf("%.3e\n", 4 * 3.14159 * (r_prime * r_prime) * deltar * DensityUnits
                                       // * SphereDensity[sphere] * pow((r_prime/(SphereCoreRadius[sphere]*LengthUnits)), -2.2));
            // }
          }
          // printf("BaryonMassOfBigSphere: %.3e Msun\n", BaryonMassOfBigSphere/SolarMass);
        }
      }
      SphereMass += BubbleTotalMass;
      printf("mass = %"GSYM", lunit = %"GSYM", dunit = %"GSYM", rho = %"GSYM", r = %"GSYM"\n",
        SphereMass, LengthUnits, DensityUnits, SphereDensity[sphere], SphereRadius[sphere]);
      break;
    
    } // ENDSWITCH SphereType
  
    // printf("\nSphere Mass (M_sun): %"FSYM"\n", SphereMass/SolarMass);
    VelocityKep = sqrt(GravConst*SphereMass/(SphereRadius[sphere]*(LengthUnits)));
    printf("Keplarian velocity: %"GSYM"\n", VelocityKep);

    if (SphereFracKeplerianRot[sphere] > 0) {
      SphereRotationalPeriod[sphere] = 2*pi*SphereRadius[sphere]*(LengthUnits)/ (SphereFracKeplerianRot[sphere]*VelocityKep);
      SphereRotationalPeriod[sphere] /= TimeUnits;
    } else
    SphereRotationalPeriod[sphere] = 0.0;

    printf("\nKeplerian Rotation Period (s): %"GSYM"\n", SphereRotationalPeriod[sphere] 
            * SphereFracKeplerianRot[sphere]*TimeUnits);
    printf("\nSphere Rotation Period (s): %"GSYM"\n", SphereRotationalPeriod[sphere]
            * TimeUnits);

    // Calculate speed of sound for this sphere
    VelocitySound[sphere] = sqrt((SphereTemperature[sphere] * Gamma * kboltz) / (mu * mh)) / VelocityUnits;
    printf("\nVelocitySound (cm s^-1): %"GSYM"\n", VelocitySound[sphere] * VelocityUnits);

  } // ENDFOR sphere

  int metal_success = 0, metal_fail = 0;

  /* Loop over grid */
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {
      /* Compute position */
      x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
      if (GridRank > 1)
      y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
      if (GridRank > 2)
      z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

      /* Loop over spheres. */
      density = InitialDensity;
      temperature = temp1 = InitialTemperature;
      sigma = sigma1 = 0;
      colour = SphereColour * CoolData.SolarMetalFractionByMass;
      metallicity = CollapseTestInitialMetallicity * CoolData.SolarMetalFractionByMass;
      HII_Fraction = CollapseTestInitialFractionHII;
      HeII_Fraction = CollapseTestInitialFractionHeII;
      HeIII_Fraction = CollapseTestInitialFractionHeIII;
      H2I_Fraction = CollapseTestInitialFractionH2I;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
        Velocity[dim] = 0;
        DMVelocity[dim] = 0;
      }

      temperature_to_intialize = InitialTemperature;
      density_to_initialize = InitialDensity;
      metallicity_to_initialize = 1.e-20;

      for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
      /* Compute Cartesian coordinates for rotational properties */  
      outer_radius = (SphereSmoothSurface[sphere] == TRUE) ? 
          SphereSmoothRadius[sphere]*SphereRadius[sphere] : SphereRadius[sphere];
      
      /* Find distance from center. */
      FLOAT xpos, ypos, zpos;
      xpos = x-SpherePosition[sphere][0];
      ypos = y-SpherePosition[sphere][1];
      zpos = z-SpherePosition[sphere][2];

      r = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
      rcyl = sqrt(xpos*xpos + ypos*ypos);
      r = max(r, 0.1*CellWidth[0][0]);
      rcyl = max(rcyl, 0.1*CellWidth[0][0]);

      if (r < outer_radius) {

        /* Compute spherical coordinate theta */
        if (xpos != 0) {
          if (xpos > 0 && ypos >= 0)      theta = atan(ypos/xpos);
          else if (xpos < 0 && ypos >= 0) theta = pi + atan(ypos/xpos);
          else if (xpos < 0 && ypos < 0)  theta = pi + atan(ypos/xpos);
          else if (xpos > 0 && ypos < 0)  theta = 2*pi + atan(ypos/xpos);
        } 
        else if (xpos == 0 && ypos > 0) theta = 3.14159 / 2.0;
        else if (xpos == 0 && ypos < 0) theta = (3*3.14159) / 2.0;
        else                            theta = 0.0;

        // Find out which shell the cell is in
        a = Ang(SphereAng1[sphere],SphereAng2[sphere],SphereRadius[sphere],r);

        /* Start with solid body rotation and then add in a
        velocity of a fraction of sound speed in a random
        direction to create turbulence */     
        if (SphereRotationalPeriod[sphere] > 0) {
          RotVelocity[0] = -2*pi*(ypos*cos(a) + zpos*sin(a)) / SphereRotationalPeriod[sphere];
          RotVelocity[1] = 2*pi*(xpos*cos(a)) / SphereRotationalPeriod[sphere];
          RotVelocity[2] = 2*pi*(xpos*sin(a)) / SphereRotationalPeriod[sphere];
          if (SphereType[sphere] == 7 || SphereType[sphere] == 9)
            GasRotVelocityCorrection = pow(rcyl/SphereRadius[sphere],-1.5);
          for (dim = 0; dim < MAX_DIMENSION; dim++) {
            Velocity[dim] = GasRotVelocityCorrection * RotVelocity[dim];
            DMVelocity[dim] = DMRotVelocityCorrection * RotVelocity[dim];
          }
        } else {
          RotVelocity[0] = RotVelocity[1] = RotVelocity[2] = 0;
        }

        Velocity[0] += SphereTurbulence[sphere] * Maxwellian(VelocitySound[sphere], VelocityUnits, mu, Gamma);
        Velocity[1] += SphereTurbulence[sphere] * Maxwellian(VelocitySound[sphere], VelocityUnits, mu, Gamma);
        Velocity[2] += SphereTurbulence[sphere] * Maxwellian(VelocitySound[sphere], VelocityUnits, mu, Gamma);
        m = 0;

        /* If collapsing uniform sphere, add radial velocity (see
        Bertschinger 1985) */
        if (SphereType[sphere] == 8) {
          ibin = (int) floor(r * NR);
          radial_velocity = vr[ibin] + 
          (vr[ibin+1] - vr[ibin]) / (radius_vr[ibin+1] - radius_vr[ibin]) *
          (r - radius_vr[ibin]);
          Velocity[0] += radial_velocity * xpos / r;
          Velocity[1] += radial_velocity * ypos / r;
          Velocity[2] += radial_velocity * zpos / r;
          DMVelocity[0] += radial_velocity * xpos / r;
          DMVelocity[1] += radial_velocity * ypos / r;
          DMVelocity[2] += radial_velocity * zpos / r;
        }
    


        /* 1) Uniform */
        if (SphereType[sphere] == 1) dens1 = SphereDensity[sphere];

        /* 2) r^-2 power law */
        if (SphereType[sphere] == 2) dens1 = SphereDensity[sphere]*pow(r/SphereRadius[sphere], -2);

        /* 3) NFW profile (use look-up table for temperature and velocity dispersion)*/
        if (SphereType[sphere] == 3) {
          x1 = r/SphereCoreRadius[sphere];
          dens1 = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));
          for (m = 1; m < NFW_POINTS; m++)
            if (r > NFWRadius[m]) {
              temp1 = NFWTemp[m] + (NFWTemp[m-1] - NFWTemp[m])*
              (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
              sigma1 = NFWSigma[m] + (NFWSigma[m-1] - NFWSigma[m])*
              (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
              break;
            }
        }

        /* 4) Gaussian */
        if (SphereType[sphere] == 4) {
          dens1 = SphereDensity[sphere] * PEXP(-0.5*pow(r/SphereCoreRadius[sphere], 2));
        }

        /* 5) r^-2 power law with core radius */
        if (SphereType[sphere] == 5) {
          if (r < SphereCoreRadius[sphere]) {
            dens1 = SphereDensity[sphere]*pow(SphereCoreRadius[sphere] / SphereRadius[sphere], -2);
          } else {
            dens1 = SphereDensity[sphere]*pow(r/SphereRadius[sphere], -2);
          }
        }

        /* 6) Bonnor-Ebert Sphere */
        if (SphereType[sphere] == 6) {
          dens1 = SphereDensity[sphere] * BE(r*Scale_Factor[sphere]);
        }

        /* 7) Uniform density, Keplerian disk */
        if (SphereType[sphere] == 7) {
          if (fabs(zpos/r) < sin_deltaDisk[sphere] && r > SphereCoreRadius[sphere]) {
            dens1 = SphereDensity[sphere];
            temperature = SphereTemperature[sphere] - GravConst * SphereMass * (mh/kboltz) * 
                (1.0/rcyl - 1.0/sqrt(rcyl*rcyl + zpos*zpos)) / LengthUnits;
            //printf("r=%"FSYM", z=%"FSYM", T=%"GSYM"\n", r,zpos,temperature);
          } else dens1 = InitialDensity;
        }

        /* 8) Uniform density, collapsing sphere after
        turnaround.  Factor of 0.79 compensates for the
        larger dispersion associated with the top-hat */
        if (SphereType[sphere] == 8) {
          dens1 = SphereDensity[sphere];
          sigma1 = SphereDispersion[sphere] * VelocityKep * r / SphereRadius[sphere];
        }

        /* 9) Circumbinary BH accretion disk (Lippai et al. 2008) */
        if (SphereType[sphere] == 9) {
          MidplaneDensity = InnerDensity * pow(rcyl/CavityRadius, -0.6);
          MidplaneTemperature = InnerTemperature * pow(rcyl/CavityRadius, -0.9);
          if (rcyl < ThickenTransitionRadius) ScaleHeight = InnerScaleHeight;
          else // (rcyl > ThickenTransitionRadius)
            ScaleHeight = InnerScaleHeight * pow(rcyl/ThickenTransitionRadius, 1.05);
          ScaleHeight = max(ScaleHeight, 2*CellWidth[0][0]);

          // printf("r=%"FSYM", z=%"FSYM", h=%"GSYM", rho=%"GSYM", T=%"GSYM"\n", 
          //      rcyl,zpos,ScaleHeight,MidplaneDensity,MidplaneTemperature);
          if (fabs(zpos) < ScaleHeight && rcyl > SphereCoreRadius[sphere] && rcyl > CavityRadius) {
            dens1 = MidplaneDensity;
            temperature = MidplaneTemperature;
          //        GravConst * SphereMass * (mh/kboltz) * 
          //        (1.0/rcyl - 1.0/sqrt(rcyl*rcyl + zpos*zpos)) / LengthUnits;
          //      printf("r=%"FSYM", z=%"FSYM", h=%"GSYM", rho=%"GSYM", T=%"GSYM"\n", 
          //       rcyl, zpos, ScaleHeight, dens1, temperature);
          } else {
            dens1 = InitialDensity;
            Velocity[0] = Velocity[1] = Velocity[2] = 0.0f;
          }
        } // ENDIF type 9
    
        /* 10) disk (ok, it's not a sphere, ... ) */
        if (SphereType[sphere] == 10) {
        FLOAT xpos, ypos, zpos, xpos1, ypos1, zpos1, zheight, drad;
        FLOAT ScaleHeightz = SphereCoreRadius[sphere]/6.0, ScaleHeightR = SphereCoreRadius[sphere];

        /* Loop over dims if using Zeus (since vel's face-centered). */

        for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0); dim++) {
          /* Compute position. */
          xpos = x-SpherePosition[sphere][0] - (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
          ypos = y-SpherePosition[sphere][1] - (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
          zpos = z-SpherePosition[sphere][2] - (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

          /* Compute z and r_perp (SphereVelocity is angular momentum 
          and must have unit length). */    
          zheight = SphereVelocity[sphere][0]*xpos + SphereVelocity[sphere][1]*ypos + SphereVelocity[sphere][2]*zpos;
          xpos1 = xpos - zheight*SphereVelocity[sphere][0];
          ypos1 = ypos - zheight*SphereVelocity[sphere][1];
          zpos1 = zpos - zheight*SphereVelocity[sphere][2];
          drad = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);

          /* If we're above the disk, then exit. */
          //    if (zheight > max(5.0*ScaleHeightz, 2.0*CellWidth[0][0]))
          //      continue;

          /* Compute density (van der Kruit & Searle 1982  1982A&A...110...61V ). */
          if (dim == 0)
          dens1 = SphereDensity[sphere]*PEXP(-drad/ScaleHeightR)/
            pow(cosh(zheight/max(ScaleHeightz, CellWidth[0][0])), 2);
          //    if (dens1 < density)
          //      break;

          /* Compute velocity magnitude (divided by drad). 
          This assumes PointSourceGravityPosition and Sphere center 
          are the same.  This should be fixed to use the disk mass
          as well, but that's a bit tricky. */
          //    float vel = sqrt(PointSourceGravityConstant/drad)/drad;
          float accel = PointSourceGravityConstant;
          if (PointSourceGravity == 1) 
            accel = PointSourceGravityConstant / (pow(drad,3) + pow(PointSourceGravityCoreRadius, 3));
          if (PointSourceGravity == 2) {
            x1 = drad/PointSourceGravityCoreRadius;
            accel = PointSourceGravityConstant*(log(1+x1)-x1/(1+x1)) / pow(drad, 3);
          }

          /* Compute velocty: L x r_perp. */
          float vel = sqrt(accel);
          if (dim == 0 || dim == 1) 
            DiskVelocity[0] = vel*(SphereVelocity[sphere][1]*zpos1 - SphereVelocity[sphere][2]*ypos1);
          if (dim == 0 || dim == 2)
            DiskVelocity[1] = vel*(SphereVelocity[sphere][2]*xpos1 + SphereVelocity[sphere][0]*zpos1);
          if (dim == 0 || dim == 3)
            DiskVelocity[2] = vel*(SphereVelocity[sphere][0]*ypos1 - SphereVelocity[sphere][1]*xpos1);
        } // end: loop over dims
        } // end: disk

        
        /************************************************************************/
        /* 11) Density and Velcoity (Bubbles Inside Sphere) */
        /* Added 2018.03.23 */

        int n_overlap = 0; // count this grid is in how many bubble regions
        if (SphereType[sphere] == 11) {
          // dens1 = SphereDensity[sphere]; // Still inside BigSphere, so set to SphereDensity[sphere] first

          if (SphereSubType[sphere] == 1) {
            if (r < SphereCoreRadius[sphere]) {
              dens1 = SphereDensity[sphere]*pow(SphereCoreRadius[sphere] / SphereRadius[sphere], SphereSubTypeSlope[sphere]);
              //fprintf(stdout,"SphereSubType[sphere]: %d", SphereSubType[sphere]);
              //fprintf(stdout,"SphereSubTypeSlope[sphere]: %f", SphereSubTypeSlope[sphere]); 
            } else {
              dens1 = SphereDensity[sphere]*pow(r/SphereRadius[sphere], SphereSubTypeSlope[sphere]);
            }
          } else if (SphereSubType[sphere] == 0) {
            dens1 = SphereDensity[sphere];
          }
                 
          float dist_to_BubbleCenter = 0.;
          FLOAT xrel, yrel, zrel;
          float overlap_density = 0.;
          
          if (MetallicityDistributionCase[sphere] == 1 || MetallicityDistributionCase[sphere] == 5) { 
            // uniform gas density for the sphere
            float bubble_radius; // 1 kpc in code units
            float inner_radius; // 200 pc in code units
            
            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              bubble_radius = BubbleRadius[sphere][count] / LengthUnits;
              inner_radius = BubbleInnerRadius[sphere][count] / LengthUnits;
              // printf("bubble_radius: %f,  inner_radius: %f\n", bubble_radius, inner_radius);

              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < inner_radius) {
                dens1 = BubbleInnerDensity[sphere][count] / DensityUnits;
              }
              else if (dist_to_BubbleCenter >= inner_radius && dist_to_BubbleCenter < bubble_radius) {
                dens1 = BubbleOuterDensity[sphere][count] / DensityUnits;
              }
              if (dist_to_BubbleCenter < bubble_radius) {
                Velocity[0] += BubbleVelocityTowardCenter[sphere][count][0];
                Velocity[1] += BubbleVelocityTowardCenter[sphere][count][1];
                Velocity[2] += BubbleVelocityTowardCenter[sphere][count][2];
              }
            }
          } else if (MetallicityDistributionCase[sphere] == 6) { 
            // uniform gas density for the sphere
            float bubble_radius; // bubble radius in code units
            float inner_radius; // bubble inner radius in code units
            
            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              bubble_radius = BubbleRadius[sphere][count] / LengthUnits;
              inner_radius = BubbleInnerRadius[sphere][count] / LengthUnits;
              // printf("bubble_radius: %f,  inner_radius: %f\n", bubble_radius, inner_radius);

              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < inner_radius) {
                dens1 = BubbleInnerDensity[sphere][count] / DensityUnits;
              }
              else if (dist_to_BubbleCenter >= inner_radius && dist_to_BubbleCenter < bubble_radius) {
                dens1 = BubbleOuterDensity[sphere][count] / DensityUnits;
              }
              if (dist_to_BubbleCenter < bubble_radius) {
                Velocity[0] += BubbleVelocityTowardCenter[sphere][count][0];
                Velocity[1] += BubbleVelocityTowardCenter[sphere][count][1];
                Velocity[2] += BubbleVelocityTowardCenter[sphere][count][2];
              }
            }
          } else if (MetallicityDistributionCase[sphere] == 2) { 
            // big bubbles
            float bubble_radius = 0.1 * outer_radius;

            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < bubble_radius) {
                Velocity[0] += BubbleVelocity[sphere][count][0];
                Velocity[1] += BubbleVelocity[sphere][count][1];
                Velocity[2] += BubbleVelocity[sphere][count][2];
                dens1 = BubbleDensity[sphere][count];
              }
            }
          } else if (MetallicityDistributionCase[sphere] == 3) { 
            // small bubbles which contain metals
            float bubble_radius = 0.05 * outer_radius;

            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < bubble_radius) {
                Velocity[0] += BubbleVelocity[sphere][count][0];
                Velocity[1] += BubbleVelocity[sphere][count][1];
                Velocity[2] += BubbleVelocity[sphere][count][2];
                dens1 = BubbleDensity[sphere][count];
              }
            }
          } else if (MetallicityDistributionCase[sphere] == 4) {
            // gradient gas density for sphere
            if (r < SphereCoreRadius[sphere]) {
              dens1 = SphereDensity[sphere];
            } else {
              dens1 = SphereDensity[sphere]*pow(r/SphereCoreRadius[sphere], -2.2);
            }

            float bubble_radius; // 1 kpc in code units
            float inner_radius; // 200 pc in code units
            
            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              bubble_radius = BubbleRadius[sphere][count] / LengthUnits;
              inner_radius = BubbleInnerRadius[sphere][count] / LengthUnits;
              // printf("bubble_radius: %f,  inner_radius: %f\n", bubble_radius, inner_radius);

              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < inner_radius) {
                dens1 = BubbleInnerDensity[sphere][count] / DensityUnits;
              }
              else if (dist_to_BubbleCenter >= inner_radius && dist_to_BubbleCenter < bubble_radius) {
                dens1 = BubbleOuterDensity[sphere][count] / DensityUnits;
              }

              if (dist_to_BubbleCenter < bubble_radius) {
                Velocity[0] += BubbleVelocityTowardCenter[sphere][count][0];
                Velocity[1] += BubbleVelocityTowardCenter[sphere][count][1];
                Velocity[2] += BubbleVelocityTowardCenter[sphere][count][2];
              }
            }
          } else if (MetallicityDistributionCase[sphere] == 7 || MetallicityDistributionCase[sphere] == 8 || MetallicityDistributionCase[sphere] == 9) { 
            // uniform gas density for the sphere
            float bubble_radius; // bubble radius in code units
            float inner_radius; // bubble inner radius in code units
            
            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              bubble_radius = BubbleRadius[sphere][count] / LengthUnits;
              inner_radius = BubbleInnerRadius[sphere][count] / LengthUnits;
              // printf("bubble_radius: %f,  inner_radius: %f\n", bubble_radius, inner_radius);

              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < inner_radius) {
                n_overlap += 1;
                overlap_density += BubbleInnerDensity[sphere][count] / DensityUnits;
              }
              else if (dist_to_BubbleCenter >= inner_radius && dist_to_BubbleCenter < bubble_radius) {
                n_overlap += 1;
                overlap_density += BubbleOuterDensity[sphere][count] / DensityUnits;
              }

              if (dist_to_BubbleCenter < bubble_radius) {
                Velocity[0] += BubbleVelocityTowardCenter[sphere][count][0];
                Velocity[1] += BubbleVelocityTowardCenter[sphere][count][1];
                Velocity[2] += BubbleVelocityTowardCenter[sphere][count][2];
              }
            }
            if (n_overlap > 3) {
              printf("n_overlap: %d\n", n_overlap);  
            }

            if (overlap_density > 0.) {
              dens1 = overlap_density;
            }

          } else if (MetallicityDistributionCase[sphere] == -1) {
            printf("WARNINR: No metal bubbles.\n");
            // return FAIL;
          } //end density for metal bubble
        } // end if SphereType[sphere] == 11 for density
        
        /* Added 2018.03.23 */
        /************************************************************************/



        /* If the density is larger than the background (or the previous
        sphere), then set the velocity. */
        if (dens1 > InitialDensity) {
        if (density <= InitialDensity) {
          old_density = 0;
          density = dens1;
        } else {
          old_density = density;
          density += dens1;
        }
        weight = dens1/density;

        if (SphereType[sphere] != 7 && SphereType[sphere] != 9) {
          if (temp1 == InitialTemperature) {
          if (SphereConstantPressure[sphere] == TRUE) {
            temperature = SphereTemperature[sphere] * (SphereDensity[sphere] / dens1);
          } else {
            temperature = SphereTemperature[sphere];
          }
          } else {
          temperature = temp1;
          }
        } // end if (SphereType[sphere] != 7 && SphereType[sphere] != 9)  

        metallicity = SphereMetallicity[sphere]; // Set to SphereMetallicity First



        /************************************************************************/
        /* 11) Bubble Metallicity and Temperature */
        /* 2018.10.18 */

        if (SphereType[sphere] == 11) {
          float dist_to_BubbleCenter = 0.;
          FLOAT xrel, yrel, zrel;                 
 
          if (SphereSubType[sphere] == 1) {
            if (r < SphereCoreRadius[sphere]) {
              temperature = 1e4;
            } else {
              temperature = SphereTemperature[sphere];
            }
          } else if (SphereSubType[sphere] == 0) {
            temperature = SphereTemperature[sphere];
          }

          if (MetallicityDistributionCase[sphere] == 0) { 
            // uniform metallicity distribution
            metallicity = 1.0e-4; // this number will multiply CoolData.SolarMetalFractionByMass                  
            temperature = temperature_to_intialize;
          } else if (MetallicityDistributionCase[sphere] == 1 || MetallicityDistributionCase[sphere] == 4 || MetallicityDistributionCase[sphere] == 5) {
            float bubble_radius; // 1 kpc radius in code units
            float inner_radius; // 200 pc radius in code units

            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              bubble_radius = BubbleRadius[sphere][count] / LengthUnits;
              inner_radius = BubbleInnerRadius[sphere][count] / LengthUnits;
              // printf("bubble_radius: %f,  inner_radius: %f\n", bubble_radius, inner_radius);

              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < inner_radius) {
              metallicity = BubbleInnerMetallicity[sphere][count];
              temperature = BubbleInnerTemperature[sphere][count];
              }
              else if (dist_to_BubbleCenter >= inner_radius && dist_to_BubbleCenter < bubble_radius) {
              float ri[9], Zi[9], Ti[9];
              for (int iter=0; iter<9; iter++) {
                ri[iter] = inner_radius + ((bubble_radius - inner_radius) / 8.) * iter;
                float zlog = 0.0;
                zlog = log10(BubbleInnerMetallicity[sphere][count]) + 
                  ((log10(1.e-6) - log10(BubbleInnerMetallicity[sphere][count]))/8.) * iter;
                Zi[iter] = pow(10, zlog); 
                // printf("Zi[%d]: %f\n", iter, Zi[iter]);
              }
              for (int iter=0; iter<8; iter++) {
                if (dist_to_BubbleCenter >= ri[iter] && dist_to_BubbleCenter < ri[iter+1]) {
                metallicity = Zi[iter];
                temperature = BubbleOuterTemperature[sphere][count];
                }
              }
              // printf("metallicity_to_initialize: %.3e\n", metallicity_to_initialize);
              // metallicity_to_initialize = BubbleOuterMetallicity[sphere][count];
              // temperature_to_intialize = BubbleOuterTemperature[sphere][count];                        
              }
            }
          } else if ( MetallicityDistributionCase[sphere] == 6) {
            float bubble_radius; // 1 kpc radius in code units
            float inner_radius; // 200 pc radius in code units

            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              bubble_radius = BubbleRadius[sphere][count] / LengthUnits;
              inner_radius = BubbleInnerRadius[sphere][count] / LengthUnits;
              // printf("bubble_radius: %f,  inner_radius: %f\n", bubble_radius, inner_radius);

              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter < inner_radius) {
              metallicity = BubbleInnerMetallicity[sphere][count];
              temperature = BubbleInnerTemperature[sphere][count];
              }
              else if (dist_to_BubbleCenter >= inner_radius && dist_to_BubbleCenter < bubble_radius) {
              float ri[9], Zi[9], Ti[9];
              for (int iter=0; iter<9; iter++) {
                ri[iter] = inner_radius + ((bubble_radius - inner_radius) / 8.) * iter;
                float zlog = 0.0;
                zlog = log10(BubbleInnerMetallicity[sphere][count]) + 
                  ((log10(1.e-6) - log10(BubbleInnerMetallicity[sphere][count]))/8.) * iter;
                Zi[iter] = pow(10, zlog); 
                // printf("Zi[%d]: %f\n", iter, Zi[iter]);
              }
              for (int iter=0; iter<8; iter++) {
                if (dist_to_BubbleCenter >= ri[iter] && dist_to_BubbleCenter < ri[iter+1]) {
                metallicity = Zi[iter];
                temperature = BubbleOuterTemperature[sphere][count];
                }
              }
              // printf("metallicity_to_initialize: %.3e\n", metallicity_to_initialize);
              // metallicity_to_initialize = BubbleOuterMetallicity[sphere][count];
              // temperature_to_intialize = BubbleOuterTemperature[sphere][count];                        
              }
            }
          } else if (MetallicityDistributionCase[sphere] == 2) { 
            // big bubbles which contain metals
            float bubble_radius = 0.1 * outer_radius;

            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter <= bubble_radius) {
              metallicity = BubbleMetallicity[sphere][count];
              temperature = BubbleTemperature[sphere][count]; //10000 K, background is 5000 K
              }
            }
          } else if (MetallicityDistributionCase[sphere] == 3) { 
            // small bubbles which contain metals
            float bubble_radius = 0.05 * outer_radius;

            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);
              
              if (dist_to_BubbleCenter <= bubble_radius) {
              metallicity = BubbleMetallicity[sphere][count];
              temperature = BubbleTemperature[sphere][count];
              }
            }
          } else if ( MetallicityDistributionCase[sphere] == 7 || MetallicityDistributionCase[sphere] == 8 || MetallicityDistributionCase[sphere] == 9 ) {
            float bubble_radius; // first in physical units then convert to in code units
            float inner_radius; // first in physical units then convert to in code units
            float Z_superposition = 0., density_superposition = 0.;

            float Z_y, Z_x, Z_a, Z_b; // function to describe metallicity profile
            float delta_Z_r, Z_r, delta_Z, linear_Z;
   
            for (int count = 0; count < NumberOfBubbles[sphere]; count++) {
              bubble_radius = BubbleRadius[sphere][count] / LengthUnits;
              inner_radius = BubbleInnerRadius[sphere][count] / LengthUnits;
              // printf("bubble_radius: %f,  inner_radius: %f\n", bubble_radius, inner_radius);

              xrel = x - center_of_bubble[sphere][count][0];
              yrel = y - center_of_bubble[sphere][count][1];
              zrel = z - center_of_bubble[sphere][count][2];
              dist_to_BubbleCenter = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);

              delta_Z_r = bubble_radius - inner_radius;
              delta_Z   = (log10(1.e-7)) - log10(BubbleInnerMetallicity[sphere][count]);
              Z_a = (delta_Z) / delta_Z_r;
              Z_b = log10(BubbleInnerMetallicity[sphere][count]);
              
              if (dist_to_BubbleCenter < inner_radius) {
                Z_superposition += BubbleInnerMetallicity[sphere][count] * BubbleInnerDensity[sphere][count] / DensityUnits;
                // metallicity_to_initialize = BubbleInnerMetallicity[sphere][count];
                temperature = BubbleInnerTemperature[sphere][count];
              }
              else if (dist_to_BubbleCenter >= inner_radius && dist_to_BubbleCenter < bubble_radius) {
                Z_y = Z_a * (dist_to_BubbleCenter - inner_radius) + Z_b;
                // printf("delta_Z_r: %.3e\n", bubble_radius - inner_radius);
                // printf("dist: %.3e\n", (dist_to_BubbleCenter - inner_radius));
                // printf("Z_a: %.3e\n", Z_a);
                // printf("Z_b: %.3e\n", Z_b);
                // printf("dist_to_BubbleCenter: %.3e\n", dist_to_BubbleCenter);
                // printf("inner_radius: %.3e\n\n\n", inner_radius);
                // printf("log10(Z_y): %.3e\n", Z_y);
                // printf("Z_y: %.3e\n", pow(10, Z_y));
                linear_Z = pow(10, Z_y);
                // printf("linear_Z: %.3e\n", linear_Z);
                // printf("BubbleOuterDensity: %.3e\n", BubbleOuterDensity[sphere][count]);

                Z_superposition += linear_Z * BubbleOuterDensity[sphere][count] / DensityUnits; 
                // printf("Z_superposition: %.3e\n", Z_superposition);
                // printf("dens1: %.3e\n", dens1);
                temperature = BubbleOuterTemperature[sphere][count];

                ////// discrete metallicity
                // float ri[9], Zi[9], Ti[9];
                // for (int iter=0; iter<9; iter++) {
                //   ri[iter] = inner_radius + ((bubble_radius - inner_radius) / 8.) * iter;
                //   float zlog = 0.0;
                //   zlog = log10(BubbleInnerMetallicity[sphere][count]) + 
                //     ((log10(1.e-6) - log10(BubbleInnerMetallicity[sphere][count]))/8.) * iter;
                //   Zi[iter] = pow(10, zlog); 
                //   // printf("Zi[%d]: %f\n", iter, Zi[iter]);
                // }
                // for (int iter=0; iter<8; iter++) {
                //   if (dist_to_BubbleCenter >= ri[iter] && dist_to_BubbleCenter < ri[iter+1]) {
                //     metallicity_to_initialize = Zi[iter];
                //     temperature_to_intialize = BubbleOuterTemperature[sphere][count];
                //   }
                // }
                // printf("metallicity_to_initialize: %.3e\n", metallicity_to_initialize);
                // metallicity_to_initialize = BubbleOuterMetallicity[sphere][count];
                // temperature_to_intialize = BubbleOuterTemperature[sphere][count]; 
                //////

              }
            } // END bubble loop
            metallicity = Z_superposition / dens1; // weighted metallicity 
            // if (metallicity_to_initialize > 0.) {
              // printf("metallicity_to_initialize: %.3e\n", metallicity_to_initialize);  
            // }
          }  else if (MetallicityDistributionCase[sphere] == -1) {
            printf("WARNING: No metal bubbles.\n");
            // return FAIL;
          } //end metallicity & temperature for metal bubble
          
          //metallicity = metallicity_to_initialize;
          
        } // end if SphereType[sphere] == 11 for density
        
        /* Added 2018.10.18 */
        /************************************************************************/
        


        sigma = sigma1;
        if (SphereType[sphere] != 10) {
          for (dim = 0; dim < GridRank; dim++)
          Velocity[dim] = Velocity[dim] + SphereVelocity[sphere][dim];
          //      Velocity[dim] = (1-weight)*Velocity[dim] + 
          //        weight*SphereVelocity[sphere][dim];
        } else // SphereType[sphere] == 10
          for (dim = 0; dim < GridRank; dim++)
          Velocity[dim] = (1-weight)*Velocity[dim] + weight*DiskVelocity[dim];

          // printf("%"ISYM" %"ISYM" %"ISYM": rho=%"GSYM", w=%"GSYM", T=%"GSYM", V=%"GSYM" %"GSYM" %"GSYM"\n", i, j, k,
          //    density, weight, temperature, Velocity[0], Velocity[1],
          //    Velocity[2]);

        /* only mark first sphere */
        if (sphere == 0)
          colour = SphereColour * CoolData.SolarMetalFractionByMass;

        HII_Fraction   = SphereHII[sphere];
        HeII_Fraction  = SphereHeII[sphere];
        HeIII_Fraction = SphereHeIII[sphere];
        H2I_Fraction   = SphereH2I[sphere];
        } // end if (dens1 > InitialDensity)

      } // end: if (r < SphereRadius)
      
      else { // if (r > SphereRadius)
        /* If collapsing uniform sphere, add radial velocity and
        density (see Bertschinger 1985) */
        if (SphereType[sphere] == 8) {
        ibin = (int) floor(r * NR);
        radial_velocity = vr[ibin] + (vr[ibin+1] - vr[ibin]) / (radius_vr[ibin+1] - 
              radius_vr[ibin]) * (r - radius_vr[ibin]);
        density = exterior_rho[ibin] + (exterior_rho[ibin+1] - exterior_rho[ibin]) / 
            (radius_vr[ibin+1] - radius_vr[ibin]) * (r - radius_vr[ibin]);

        // IGM in pressure equilibrium (2.388 is the density at
        // the turnaround radius)
        temperature = InitialTemperature / (density / 2.388);

        Velocity[0] += radial_velocity * xpos / r;
        Velocity[1] += radial_velocity * ypos / r;
        Velocity[2] += radial_velocity * zpos / r;
        DMVelocity[0] += radial_velocity * xpos / r;
        DMVelocity[1] += radial_velocity * ypos / r;
        DMVelocity[2] += radial_velocity * zpos / r;

        // Add rotational velocities
        if (SphereRotationalPeriod[sphere] > 0) {
          RotVelocity[0] = -2*pi*(ypos*cos(a) + zpos*sin(a)) / SphereRotationalPeriod[sphere];
          RotVelocity[1] = 2*pi*(xpos*cos(a)) / SphereRotationalPeriod[sphere];
          RotVelocity[2] = 2*pi*(xpos*sin(a)) / SphereRotationalPeriod[sphere];
          for (dim = 0; dim < MAX_DIMENSION; dim++) {
          DMVelocity[dim] += DMRotVelocityCorrection * RotVelocity[dim];
          Velocity[dim] += GasRotVelocityCorrection * RotVelocity[dim];
          }
        }
        } // ENDIF type == 8
      } // end: if (r > SphereRadius)

      if (SphereSmoothSurface[sphere] == TRUE && r < SphereSmoothRadius[sphere]*SphereRadius[sphere] &&
        r > SphereRadius[sphere]) {
        
        float ramp = 1.0 - 1.0 * tanh((3.0/(SphereSmoothRadius[sphere]-1.0)) * (r/SphereRadius[sphere] - 1.0));
        ramp = max(ramp, 1.0/density);
        density *= ramp;
        if (SphereConstantPressure[sphere] == TRUE) {
        temperature /= ramp;
        }
      } // end: if (SmoothSurface)
      } // end: loop over spheres

      /* Set density. */
      BaryonField[0][n] = density*BaryonMeanDensity;

      if (StarMakerTypeIaSNe) BaryonField[MetalIaNum][n] = 1.0e-20 * BaryonField[0][n];
      // Added * BaryonField[0][n] 2018.03.15

      /* If doing multi-species (HI, etc.), set these. */
      if (MultiSpecies > 0) {

      BaryonField[HIINum][n] = HII_Fraction * CoolData.HydrogenFractionByMass * BaryonField[0][n] *
            sqrt(OmegaMatterNow) / (OmegaMatterNow*BaryonMeanDensity*HubbleConstantNow);
      //      (CosmologySimulationOmegaBaryonNow*HubbleConstantNow);

      BaryonField[HeIINum][n]  = HeII_Fraction * BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
      BaryonField[HeIIINum][n] = HeIII_Fraction * BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
      BaryonField[HeINum][n]   = (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] - 
              BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];

      if (MultiSpecies > 1) {
        BaryonField[HMNum][n]   = CollapseTestInitialFractionHM * BaryonField[HIINum][n]* pow(temperature,float(0.88));
        BaryonField[H2IINum][n] = CollapseTestInitialFractionH2II * 2.0*BaryonField[HIINum][n]* pow(temperature,float(1.8));
        if (ComovingCoordinates) 
        BaryonField[H2INum][n] = H2I_Fraction * BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1)*
              pow(OmegaMatterNow, float(1.5)) / (OmegaMatterNow*BaryonMeanDensity)/
              //        CosmologySimulationOmegaBaryonNow/
              HubbleConstantNow*2.0;
        else
        BaryonField[H2INum][n] = H2I_Fraction * CoolData.HydrogenFractionByMass * BaryonField[0][n];
      }

      BaryonField[HINum][n] = CoolData.HydrogenFractionByMass*BaryonField[0][n] - BaryonField[HIINum][n];
      if (MultiSpecies > 1)
        BaryonField[HINum][n] -= BaryonField[HMNum][n] + BaryonField[H2IINum][n] + BaryonField[H2INum][n];

      BaryonField[DeNum][n] = BaryonField[HIINum][n] + 0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
      if (MultiSpecies > 1)
        BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] - BaryonField[HMNum][n];

      /* Set Deuterium species (assumed to be negligible). */
      if (MultiSpecies > 2) {
        BaryonField[DINum][n]  = CoolData.DeuteriumToHydrogenRatio * BaryonField[HINum][n];
        BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio * BaryonField[HIINum][n];
        BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio * BaryonField[H2INum][n];
      }
      }

      /* If there are metals, set it. */
      if (SphereUseMetals) {
        BaryonField[MetalNum][n] = metallicity * CoolData.SolarMetalFractionByMass * BaryonField[0][n];
        if (metallicity >= 1.e-4) {
          // printf("Metal_Density: %.3e\n", BaryonField[MetalNum][n]);
        }
      }

      /* If there is a colour field, set it. */
      if (SphereUseColour)
      BaryonField[ColourNum][n] = colour * BaryonField[0][n];

      /* Set Velocities. */
      for (dim = 0; dim < GridRank; dim++)
      BaryonField[ivel+dim][n] = Velocity[dim] + UniformVelocity[dim];

      /* Set energy (thermal and then total if necessary). */
      BaryonField[1][n] = temperature/TemperatureUnits/ ((Gamma-1.0)*mu);

      if (DualEnergyFormalism)
      BaryonField[2][n] = BaryonField[1][n];

      if (HydroMethod != Zeus_Hydro)
      for (dim = 0; dim < GridRank; dim++)
        BaryonField[1][n] += 0.5*pow(BaryonField[ivel+dim][n], 2);

      /* Set particles if being used (generate a number of particle
      proportional to density). */

      if (SphereUseParticles)
      if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
        j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
        k >= GridStartIndex[2] && k <= GridEndIndex[2]  ) {
      
        ParticleCount += density/pow(float(RefineBy), GridRank*level);
        while (ParticleCount > 1) {
        if (SetupLoopCount > 0) {
          ParticleMass[npart] = ParticleMeanDensity * pow(float(RefineBy), GridRank*level);
          ParticleNumber[npart] = CollapseTestParticleCount++;
          ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;

          /* Set random position within cell. */
          ParticlePosition[0][npart] = x + CellWidth[0][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
          ParticlePosition[1][npart] = y + CellWidth[1][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
          ParticlePosition[2][npart] = z + CellWidth[2][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);

          /* Set bulk velocity. */
          for (sphere = 0; sphere<NumberOfSpheres;sphere++) {
          xdist = ParticlePosition[0][npart]-SpherePosition[sphere][0];
          ydist = ParticlePosition[1][npart]-SpherePosition[sphere][1];
          zdist = ParticlePosition[2][npart]-SpherePosition[sphere][2];
          for (dim = 0; dim < GridRank; dim++)
            if (sqrt(xdist*xdist + ydist*ydist + zdist*zdist) < SphereRadius[sphere]) { //particle is inside sphere
            ParticleVelocity[dim][npart] = DMVelocity[dim]+UniformVelocity[dim]+SphereVelocity[sphere][dim]; 
            } else { //particle is outside sphere
            ParticleVelocity[dim][npart] = DMVelocity[dim]+UniformVelocity[dim]; 
            }
          }
          
          /* Add random velocity; */
          if (sigma != 0) 
          for (dim = 0; dim < GridRank; dim++)
            ParticleVelocity[dim][npart] += gasdev()*sigma/VelocityUnits;
        }
        npart++;
        ParticleCount -= 1.0;
        }
      } // end: if statement

  } // end loop over grid
  // printf("Metal SN  Enriched: %d\n", metal_success);
  // printf("Metal Not Enriched: %d\n", metal_fail);

  } // end loop SetupLoopCount

  if (SphereUseParticles && debug)
  printf("CollapseTestInitialize: NumberOfParticles = %"ISYM"\n", NumberOfParticles);

  return SUCCESS;
}


/* Routine to return a Gaussian random deviate with zero mean and unite
   variance (adapted from Numerical Recipes). */

static int gasdev_iset = 0;
static float gasdev_gset;

float gasdev()
{
  float v1, v2, r = 0, fac, gasdev_ret;
  if (gasdev_iset == 0) {
  while (r >= 1 || r == 0) {
    v1 = 2.0*float(rand())/(float(RAND_MAX)) - 1.0;
    v2 = 2.0*float(rand())/(float(RAND_MAX)) - 1.0;
    r = v1*v1 + v2*v2;
  }
  fac = sqrt(-2.0*log(r)/r);
  gasdev_gset = v1*fac;
  gasdev_ret  = v2*fac;
  gasdev_iset = 1;
  } else {
  gasdev_ret  = gasdev_gset;
  gasdev_iset = 0;
  }
  return gasdev_ret;
}

/************************************************************************/

double BE(double r)
{
  double factor;
  factor = 4.8089e-04*pow(r,5) - 1.0173e-02*pow(r,4) + 7.7899e-02*pow(r,3) - 
    2.3299e-01*pow(r,2) + 1.4721e-02*r + 1.0008e+00;
  return factor;
}

/************************************************************************/

double q(double r)
{
  double factor;
  factor = 0.0015970*pow(r,5) - 0.0229113*pow(r,4) + 0.0386709*pow(r,3) + 
    0.7350457*pow(r,2) - 0.5490283*r + 0.0872061;
  return factor;
}

/************************************************************************/

double Ang(double a1, double a2, double R, double r)
{
  return ((a2-a1)/R)*r + a1;
}

/************************************************************************/

double Maxwellian(double c_tilda, double vel_unit, double mu, double gamma)
{
  // Set constants
  double mh = 1.67e-24;
  double kboltz = 1.38e-16;
  double pi = 3.14159;

  // Compute temperature in cgs units
  double c = c_tilda*vel_unit;
  double T = (pow(c,2)*mu*mh)/(kboltz*gamma);

  // Compute random Maxwellian velocity
  double mean = 0;
  double stdev = sqrt((kboltz*T)/(mu*mh));
  double u1 = rand();
  u1 = u1/RAND_MAX;
  double u2 = rand();
  u2 = u2/RAND_MAX;
  double x1 = mean + stdev*sqrt(-2*log(u1))*cos(2*pi*u2);
  double x2 = mean + stdev*sqrt(-2*log(u1))*sin(2*pi*u2);
  return (x1/vel_unit);
}

double ERF(double x)
{
  return (2.0 / sqrt(M_PI)) * (x - pow(x,3)/3.0 + pow(x,5)/10.0 - pow(x,7)/42.0 + pow(x,9)/216.0);
}

int ComputeRadialVelocity(float density, double mass, float r_init, 
        double VelocityUnits, double LengthUnits,
        float dm_rotvel_corr, double radius_vr[], double Vr[], 
        double exterior_rho[], int Npts)
{

  const float theta0 = 0.5*M_PI, theta1 = 1.9*M_PI;
  const float Grav = 6.673e-8, Mpc = 3.086e24, yr = 3.1557e7, Msun = 1.989e33;
  const float delta_i = 0.5;
  float dtheta, delta_t[NTHETA], Theta[NTHETA], Lambda[NTHETA], Chi;
  float beta[NTHETA], little_d[NTHETA];
  float local_Vr[NTHETA], local_radius[NTHETA], local_rho[NTHETA], vmax, Vl;
  float t_init, t_i, r_i, dr;
  float z_vir, t_vir;
  float factor, t_ta, rho_ci, r_ta_now;
  int i, n, iinit, iturn, i_boundary, ithis;

  dtheta = (theta1 - theta0) / (NTHETA-1.0);
  for (i = 0; i < NTHETA; i++) {
    Theta[i] = theta0 + i*dtheta;
    Lambda[i] = pow(sin(Theta[i]/2), 2) * pow((Theta[i] - sin(Theta[i])) / M_PI, -8.0/9);
    delta_t[i] = 4.5 * pow((Theta[i] - sin(Theta[i])), 2) / pow((1.0 - cos(Theta[i])), 3) - 1.0;
    beta[i] = sin(0.5*Theta[i]) * sin(0.5*Theta[i]);
    little_d[i] = 0.75 * (Theta[i] - sin(Theta[i]));
    if (delta_t[i] < density)
      iinit = i;
    if (Theta[i] < M_PI)   // Turnaround
      iturn = i;
  }

  t_init = 5.38e8 * yr * pow((1+InitialRedshift) / 10.0, -1.5);
  z_vir  = pow((Theta[iinit] - sin(Theta[iinit])) / (2*M_PI), 2.0/3) * (1 + InitialRedshift) - 1.0;
  t_vir  = 5.38e8 * yr * pow((1+z_vir) / 10.0, -1.5);
  t_ta   = 0.5 * t_vir;
  t_i    = t_ta * pow(delta_i, 1.5) / (3*M_PI/4);
  rho_ci = 1.0 / (6 * M_PI * Grav * t_i*t_i);
  r_i    = pow(mass / (4*M_PI/3 * rho_ci), 1.0/3);

  dm_rotvel_corr = 2 * pow(2*r_init/r_i, 2);
  //  dm_rotvel_corr = 1;

  r_ta_now = pow((3*M_PI/4), -8.0/9) * pow(delta_i, 1.0/3) * r_i * pow(t_init / t_i, 8.0/9);

  for (i = 0; i < NTHETA; i++) {
    local_Vr[NTHETA-i-1] = (Lambda[i] * r_ta_now / t_init) * sin(Theta[i]) * 
          (Theta[i] - sin(Theta[i])) / pow(1.0 - cos(Theta[i]), 2);
    Vl = local_Vr[NTHETA-i-1] / (r_ta_now / t_init);

    // Subtract Hubble flow and convert to code units
    local_Vr[NTHETA-i-1]     = (Vl - (2.0/3)*Lambda[i]) * r_ta_now / t_init;
    local_Vr[NTHETA-i-1]    /= VelocityUnits;
    Chi                      = 1.0 - 1.5 * Vl / Lambda[i];
    local_rho[NTHETA-i-1]    = pow(little_d[i],2) / pow(beta[i],3) / (1.0 + 3*Chi);
    local_radius[NTHETA-i-1] = Lambda[i] * r_ta_now / LengthUnits;
    if (local_radius[NTHETA-i-1] > r_init)
      i_boundary = NTHETA-i-1;
  }

  dr = 1.0 / Npts;
  vmax = local_Vr[i_boundary];
  for (i = 0; i < Npts; i++) {
    radius_vr[i] = i*dr;
    if (radius_vr[i] > r_init) {
      // Find bin (slow, but we only do this in the initialization, so no worries.
      for (n = 0; n < NTHETA-1; n++)
        if ((radius_vr[i] >  local_radius[n] && radius_vr[i] <= local_radius[n+1]) || n == NTHETA-2) {
          ithis = n;
          break;
        }

      // Interpolate
      Vr[i] = local_Vr[ithis] + (local_Vr[ithis+1] - local_Vr[ithis]) /
        (local_radius[ithis+1] - local_radius[ithis]) * (radius_vr[i] - local_radius[ithis]);
      exterior_rho[i] = local_rho[ithis] + (local_rho[ithis+1] - local_rho[ithis]) /
          (local_radius[ithis+1] - local_radius[ithis]) * (radius_vr[i] - local_radius[ithis]);
    } else {
      Vr[i] = radius_vr[i] / r_init * vmax;
      exterior_rho[i] = density;
    }
    //    printf("%"ISYM": r = %"GSYM" pc, v_r = %"GSYM" km/s, rho = %"GSYM"\n", 
    //     i, radius_vr[i]*LengthUnits/3.086e18, Vr[i]*VelocityUnits/1e5,
    //     exterior_rho[i]);
  }
  return SUCCESS;
}

/* function of velocity toward center */
double ComputeBubbleVelocityInducedByDarkMatter(const float *PointSourceGravityPosition, 
    const float PointSourceGravityCoreRadius, const float *center_of_bubble, 
    float *BubbleVelocityTowardCenter, const float DensityUnits, const float LengthUnits, 
    const float TemperatureUnits, const float TimeUnits, const float VelocityUnits, const float SolarMass,
    const float GravConst, const float SphereDensity, const int MetallicityDistributionCase, const float SphereCoreRadius) {
  
  // printf("MetallicityDistributionCase: %d\n", MetallicityDistributionCase);
  float pc_in_cm = 3.0857e+18, steps = 1.e6, rho_0, G_physical = 6.67e-8; // cgs
  int iter_steps = 1e6;
  // PointSourceGravityConstant: 1e8 Msun...
  rho_0 = PointSourceGravityConstant / (4 * 3.14159 * pow(PointSourceGravityCoreRadius,3) / 3 );
  // printf("rho_0: %.3e\n", rho_0);

  // printf("coreradius: %.3e\n", PointSourceGravityCoreRadius);
  // printf("coremass: %.3e\n", PointSourceGravityConstant);

  float dist_to_halocenter, xx, yy, zz, norm_factor;
  xx = center_of_bubble[0] - PointSourceGravityPosition[0]; // code unit
  // printf("center_of_bubble[0]: %f, PointSourceGravityPosition[0]: %f\n", center_of_bubble[0], PointSourceGravityPosition[0]);
  yy = center_of_bubble[1] - PointSourceGravityPosition[1]; // code unit
  zz = center_of_bubble[2] - PointSourceGravityPosition[2]; // code unit
  norm_factor = sqrt(xx*xx + yy*yy + zz*zz);
  dist_to_halocenter = sqrt(xx*xx + yy*yy + zz*zz);
  // printf("dist_to_halocenter: %.3e code unit\n", dist_to_halocenter);
  // printf("dist_to_halocenter: %.3e pc\n", dist_to_halocenter*LengthUnits/pc_in_cm);
  
  float core_ratio, vmag, DMmass_within_radius;
  core_ratio = PointSourceGravityCoreRadius / LengthUnits; // core radius in code unit

  if (dist_to_halocenter <= core_ratio) {
    DMmass_within_radius = PointSourceGravityConstant * pow((dist_to_halocenter/core_ratio), 3); // cgs
    vmag = sqrt(G_physical * DMmass_within_radius / dist_to_halocenter);
    vmag /= VelocityUnits;
  }
  else if (dist_to_halocenter > core_ratio) {
    DMmass_within_radius = PointSourceGravityConstant; // cgs
    // printf("1e6 M_sun: %.3e M_sun\n", DMmass_within_radius/SolarMass);
    
    float r_prime = PointSourceGravityCoreRadius; // code unit
    // printf("coreradius: %.3e\n", r_prime);
    
    float deltar = (dist_to_halocenter*LengthUnits - PointSourceGravityCoreRadius) / 1000000.; // physcial length
    // printf("deltar: %.3e\n", deltar);
    
    for (int iteration = 0; iteration < 1000000; iteration ++ ) {
      r_prime += deltar;
      DMmass_within_radius += 4 * 3.14159 * (r_prime * r_prime) * deltar * rho_0 * pow((r_prime/PointSourceGravityCoreRadius), -2.2);

      if (iteration % 100000 == 0) {
        // printf("r_prime: %.3e pc\n", r_prime/pc_in_cm);
        // printf("DarkMatterMass: %.3e M_sun\n", DMmass_within_radius/SolarMass);
      }      
    }
  }
  // printf("Dark Matter Mass: %.3e M_sun\n", DMmass_within_radius/SolarMass);

  float BaryonMassWithinRadius, dens1, SphereCoreBaryonMass;

  SphereCoreBaryonMass = 4. * 3.14159 * pow(SphereCoreRadius*LengthUnits, 3) / 3. * SphereDensity * DensityUnits;
  // printf("SphereCoreBaryonMass: %.3e Msun\n", SphereCoreBaryonMass/SolarMass);
  // printf("SphereCoreRadius: %.3e\n", SphereCoreRadius);

  if (MetallicityDistributionCase == 1) {
    BaryonMassWithinRadius = 4. * 3.14159 * pow(dist_to_halocenter*LengthUnits, 3) / 3. * SphereDensity * DensityUnits;
    // printf("BaryonMassWithinRadius: %.3e Msun\n", BaryonMassWithinRadius/SolarMass);
  }
  else if (MetallicityDistributionCase >= 4) {
    if (dist_to_halocenter < SphereCoreRadius) { // already in code unit
      BaryonMassWithinRadius = 4. * 3.14159 * pow(dist_to_halocenter*LengthUnits, 3) / 3. * SphereDensity * DensityUnits;
    } else if (dist_to_halocenter >= SphereCoreRadius) {
      BaryonMassWithinRadius = SphereCoreBaryonMass;
      
      float r_prime = SphereCoreRadius*LengthUnits; // code unit to cgs
      float deltar = (dist_to_halocenter*LengthUnits - SphereCoreRadius*LengthUnits) / 1000000.;
      
      for (int iteration = 0; iteration < 1000000; iteration ++ ) {
        r_prime += deltar;
        BaryonMassWithinRadius += 4 * 3.14159 * (r_prime * r_prime) * deltar * DensityUnits
                                   * SphereDensity * pow((r_prime/(SphereCoreRadius*LengthUnits)), -2.2);
        if (iteration % 100000 == 0) {
          // printf("r_prime: %.3e pc\n", r_prime/pc_in_cm);
          // printf("BaryonMass: %.6e M_sun\n", BaryonMassWithinRadius/SolarMass);
        }
      }
    }
    // printf("BaryonMassWithinRadius: %.3e Msun\n", BaryonMassWithinRadius/SolarMass);
  }

  float TotalMass;
  TotalMass = DMmass_within_radius + BaryonMassWithinRadius;
  if (MyProcessorNumber == ROOT_PROCESSOR) {
     printf("TotalMass: %.3e Msun\n\n\n\n", TotalMass/SolarMass);
     printf("distance: %.3e\n", dist_to_halocenter*LengthUnits);
     printf("G: %.3e\n", G_physical);
   }

  vmag = sqrt(G_physical * TotalMass / (dist_to_halocenter*LengthUnits));
  //printf("vmag: %.3e km/s\n", vmag/100000.);    
  vmag /= VelocityUnits;
  // printf("vmag: %.3e code unit\n\n\n\n\n\n\n", vmag);

  return vmag;
}
