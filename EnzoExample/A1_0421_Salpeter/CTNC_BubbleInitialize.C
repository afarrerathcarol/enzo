/***********************************************************************
/
/  INITIALIZE A COLLAPSE TEST
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:
/
/  PURPOSE:
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static float CollapseTestInitialFractionHII   = 1.2e-5;
static float CollapseTestInitialFractionHeII  = 1.0e-14;
static float CollapseTestInitialFractionHeIII = 1.0e-17;
static float CollapseTestInitialFractionHM    = 2.0e-9;
static float CollapseTestInitialFractionH2I   = 2.0e-20;
static float CollapseTestInitialFractionH2II  = 3.0e-14;

int main()
{

  FILE *fptr;
  fptr = fopen("CollapseTestNonCosmological.enzo", "r");
  if (fptr) printf("Successfully read in CollapseTestNonCosmological.enzo.\n");
  else printf("Fail read in CollapseTestNonCosmological.enzo.\n");
  
  printf("DON'T FORGET TO CHECK FOR THE MAX_***.\n");

  /* declarations */

  int MAX_SPHERES = 10, MAX_DIMENSION = 3, MAX_BUBBLES = 100, MAX_LINE_LENGTH = 10000;
  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i, bubble;

  /* set default parameters */

  int CollapseTestNumberOfSpheres = 1;
  int CollapseTestRefineAtStart   = 1;
  int CollapseTestUseParticles    = 0;
  float CollapseTestParticleMeanDensity = 0.0;
  int CollapseTestUseColour       = 0;
  int CollapseTestUseMetals       = 0;
  float CollapseTestInitialTemperature = 1000;
  float CollapseTestInitialDensity     = 1.0;
  float CollapseTestSphereDensity[MAX_SPHERES],
  CollapseTestSphereTemperature[MAX_SPHERES],
  CollapseTestSphereVelocity[MAX_SPHERES][MAX_DIMENSION],
  CollapseTestUniformVelocity[MAX_DIMENSION],
  CollapseTestFracKeplerianRot[MAX_SPHERES],
  CollapseTestSphereTurbulence[MAX_SPHERES],
  CollapseTestSphereDispersion[MAX_SPHERES],
  CollapseTestSphereCutOff[MAX_SPHERES],
  CollapseTestSphereAng1[MAX_SPHERES],
  CollapseTestSphereAng2[MAX_SPHERES],
  CollapseTestSphereMetallicity[MAX_SPHERES],
  CollapseTestSphereSmoothRadius[MAX_SPHERES],
  CollapseTestSphereHIIFraction[MAX_SPHERES],
  CollapseTestSphereHeIIFraction[MAX_SPHERES],
  CollapseTestSphereHeIIIFraction[MAX_SPHERES],
  CollapseTestSphereH2IFraction[MAX_SPHERES],
  /* 2018.02.01 added */
  CollapseTestSphereHMFraction[MAX_SPHERES],
  CollapseTestSphereH2IIFraction[MAX_SPHERES];
  /* 2018.02.01 added */
  int CollapseTestSphereNumShells[MAX_SPHERES],
  CollapseTestSphereInitialLevel[MAX_SPHERES],
  CollapseTestSphereType[MAX_SPHERES],
  CollapseTestSphereConstantPressure[MAX_SPHERES],
  CollapseTestSphereSmoothSurface[MAX_SPHERES];
  float CollapseTestSphereRadius[MAX_SPHERES],
  CollapseTestSphereCoreRadius[MAX_SPHERES],
  CollapseTestSpherePosition[MAX_SPHERES][MAX_DIMENSION];
  float LengthUnits;

  /* 2018.03.18 added */
  int NumberOfBubbles[MAX_SPHERES], 
    MetallicityDistributionCase[MAX_SPHERES];
  /* 2018.03.18 added */  

  /* 2018.03.23 */
  float center_of_bubble[MAX_SPHERES][MAX_BUBBLES][MAX_DIMENSION];
  float BubbleVelocity[MAX_SPHERES][MAX_BUBBLES][MAX_DIMENSION], BubbleDensity[MAX_SPHERES][MAX_BUBBLES];
  float BubbleMetallicity[MAX_SPHERES][MAX_BUBBLES], BubbleTemperature[MAX_SPHERES][MAX_BUBBLES];
  // For MetallicityDistributionCase = 1
  float BubbleInnerDensity[MAX_SPHERES][MAX_BUBBLES], BubbleOuterDensity[MAX_SPHERES][MAX_BUBBLES];
  float BubbleInnerMetallicity[MAX_SPHERES][MAX_BUBBLES], BubbleOuterMetallicity[MAX_SPHERES][MAX_BUBBLES];
  float BubbleInnerTemperature[MAX_SPHERES][MAX_BUBBLES], BubbleOuterTemperature[MAX_SPHERES][MAX_BUBBLES];
  float BubbleRadius[MAX_SPHERES][MAX_BUBBLES], BubbleInnerRadius[MAX_SPHERES][MAX_BUBBLES];
  /* 2018.03.23 */

  int CollapseTestSphereVelocityTowardCenter[MAX_SPHERES], CollapseTestSphereDarkMatterHaloMass=0; 

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    CollapseTestSphereRadius[sphere]     = 1.0;
    CollapseTestSphereCoreRadius[sphere] = 0.1;
    CollapseTestSphereDensity[sphere]    = 1.0;
    CollapseTestSphereTemperature[sphere] = 1.0;
    CollapseTestFracKeplerianRot[sphere] = 0.0;
    CollapseTestSphereTurbulence[sphere] = 0.0;
    CollapseTestSphereDispersion[sphere] = 0.0;
    CollapseTestSphereCutOff[sphere] = 6.5;
    CollapseTestSphereAng1[sphere] = 0;
    CollapseTestSphereAng2[sphere] = 0;
    CollapseTestSphereNumShells[sphere] = 1;
    CollapseTestSphereSmoothRadius[sphere] = 1.2;
    CollapseTestSphereMetallicity[sphere] = 1.0e-20;
    CollapseTestSphereInitialLevel[sphere] = 0;
    CollapseTestSphereHIIFraction[sphere] = CollapseTestInitialFractionHII;
    CollapseTestSphereHeIIFraction[sphere] = CollapseTestInitialFractionHeII;
    CollapseTestSphereHeIIIFraction[sphere] = CollapseTestInitialFractionHeIII;
    CollapseTestSphereH2IFraction[sphere] = CollapseTestInitialFractionH2I;
    /* 2018.02.01 added */
    CollapseTestSphereHMFraction[sphere] = CollapseTestInitialFractionHM;
    CollapseTestSphereH2IIFraction[sphere] = CollapseTestInitialFractionH2II;
    /* 2018.02.01 added */

    /* 2018.03.18 added */
    NumberOfBubbles[sphere] = 0;
    MetallicityDistributionCase[sphere] = -1;
    /* 2018.03.18 added */

    CollapseTestSphereVelocityTowardCenter[sphere] = 0;

    for (int ii=0; ii < MAX_BUBBLES; ii++) {
      BubbleDensity[sphere][ii] = 1.0;
      BubbleInnerDensity[sphere][ii] = 1.0;
      BubbleOuterDensity[sphere][ii] = 1.0;

      BubbleInnerTemperature[sphere][ii] = 1.e4;
      BubbleOuterTemperature[sphere][ii] = 1.e4;

      BubbleInnerMetallicity[sphere][ii] = 1.e-4;
      BubbleOuterMetallicity[sphere][ii] = 1.e-4;

      BubbleRadius[sphere][ii] = 0.0;
      BubbleInnerRadius[sphere][ii] = 0.0;
    }

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CollapseTestSpherePosition[sphere][dim] = 0.5;
      CollapseTestSphereVelocity[sphere][dim] = 0;

      for (int ii =0; ii < MAX_BUBBLES; ii++) {
      center_of_bubble[sphere][ii][dim] = 0.5;
      BubbleVelocity[sphere][ii][dim] = 0.0;
      }
    }
    CollapseTestSphereType[sphere]       = 0;
    CollapseTestSphereConstantPressure[sphere] = 0;
    CollapseTestSphereSmoothSurface[sphere] = 0;
  }

  for (dim = 0; dim < MAX_DIMENSION; dim++)
  CollapseTestUniformVelocity[dim] = 0;



  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "CollapseTestNumberOfSpheres = %d",
      &CollapseTestNumberOfSpheres);
    ret += sscanf(line, "CollapseTestRefineAtStart = %d", 
      &CollapseTestRefineAtStart);
    ret += sscanf(line, "CollapseTestUseParticles = %d", 
      &CollapseTestUseParticles);
    ret += sscanf(line, "CollapseTestParticleMeanDensity = %f",
      &CollapseTestParticleMeanDensity);
    ret += sscanf(line, "CollapseTestUseColour = %d", 
      &CollapseTestUseColour);
    ret += sscanf(line, "CollapseTestUseMetals = %d", 
      &CollapseTestUseMetals);
    ret += sscanf(line, "CollapseTestInitialTemperature = %f", 
      &CollapseTestInitialTemperature);
    ret += sscanf(line, "CollapseTestInitialDensity = %f",
      &CollapseTestInitialDensity);
    ret += sscanf(line, "CollapseTestUniformVelocity = %f %f %f", 
      CollapseTestUniformVelocity, CollapseTestUniformVelocity+1,
      CollapseTestUniformVelocity+2);
    if (sscanf(line, "CollapseTestSphereType[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereType[%d] = %d", &sphere,
        &CollapseTestSphereType[sphere]);
    if (sscanf(line, "CollapseTestSphereRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereRadius[%d] = %f", &sphere,
        &CollapseTestSphereRadius[sphere]);
    if (sscanf(line, "CollapseTestSphereCoreRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereCoreRadius[%d] = %f", &sphere,
        &CollapseTestSphereCoreRadius[sphere]);
    if (sscanf(line, "CollapseTestSphereDensity[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereDensity[%d] = %f", &sphere,
        &CollapseTestSphereDensity[sphere]);
    if (sscanf(line, "CollapseTestSphereTemperature[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereTemperature[%d] = %f", &sphere,
        &CollapseTestSphereTemperature[sphere]);
    if (sscanf(line, "CollapseTestSphereMetallicity[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereMetallicity[%d] = %f", &sphere,
        &CollapseTestSphereMetallicity[sphere]);
    if (sscanf(line, "CollapseTestSpherePosition[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSpherePosition[%d] = %f %f %f", 
        &sphere, &CollapseTestSpherePosition[sphere][0],
        &CollapseTestSpherePosition[sphere][1],
        &CollapseTestSpherePosition[sphere][2]);
    if (sscanf(line, "CollapseTestSphereVelocity[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereVelocity[%d] = %f %f %f", 
        &sphere, &CollapseTestSphereVelocity[sphere][0],
        &CollapseTestSphereVelocity[sphere][1],
        &CollapseTestSphereVelocity[sphere][2]);
    if (sscanf(line, "CollapseTestFracKeplerianRot[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestFracKeplerianRot[%d] = %f", &sphere,
          &CollapseTestFracKeplerianRot[sphere]);
    if (sscanf(line, "CollapseTestSphereTurbulence[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereTurbulence[%d] = %f", &sphere,
          &CollapseTestSphereTurbulence[sphere]);
    if (sscanf(line, "CollapseTestSphereDispersion[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereDispersion[%d] = %f", &sphere,
          &CollapseTestSphereDispersion[sphere]);
    if (sscanf(line, "CollapseTestSphereCutOff[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereCutOff[%d] = %f", &sphere,
          &CollapseTestSphereCutOff[sphere]);
    if (sscanf(line, "CollapseTestSphereAng1[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereAng1[%d] = %f", &sphere,
          &CollapseTestSphereAng1[sphere]);
    if (sscanf(line, "CollapseTestSphereAng2[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereAng2[%d] = %f", &sphere,
          &CollapseTestSphereAng2[sphere]);
    if (sscanf(line, "CollapseTestSphereNumShells[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereNumShells[%d] = %d", &sphere,
          &CollapseTestSphereNumShells[sphere]);
    if (sscanf(line, "CollapseTestSphereInitialLevel[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereInitialLevel[%d] = %d", &sphere,
          &CollapseTestSphereInitialLevel[sphere]);
    if (sscanf(line, "CollapseTestSphereConstantPressure[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereConstantPressure[%d] = %d", &sphere,
        &CollapseTestSphereConstantPressure[sphere]);
    if (sscanf(line, "CollapseTestSphereSmoothSurface[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereSmoothSurface[%d] = %d", &sphere,
        &CollapseTestSphereSmoothSurface[sphere]);
    if (sscanf(line, "CollapseTestSphereSmoothRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereSmoothRadius[%d] = %f", &sphere,
        &CollapseTestSphereSmoothRadius[sphere]);
    if (sscanf(line, "CollapseTestSphereHIIFraction[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereHIIFraction[%d] = %f", &sphere,
          &CollapseTestSphereHIIFraction[sphere]);
    if (sscanf(line, "CollapseTestSphereHeIIFraction[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereHeIIFraction[%d] = %f", &sphere,
          &CollapseTestSphereHeIIFraction[sphere]);
    if (sscanf(line, "CollapseTestSphereHeIIIFraction[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereHeIIIFraction[%d] = %f", &sphere,
          &CollapseTestSphereHeIIIFraction[sphere]);
    if (sscanf(line, "CollapseTestSphereH2IFraction[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereH2IFraction[%d] = %f", &sphere,
          &CollapseTestSphereH2IFraction[sphere]);
    /* 2018.02.01 added */
    if (sscanf(line, "CollapseTestSphereHMFraction[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereHMFraction[%d] = %f", &sphere,
          &CollapseTestSphereHMFraction[sphere]);
    if (sscanf(line, "CollapseTestSphereH2IIFraction[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereH2IIFraction[%d] = %f", &sphere,
          &CollapseTestSphereH2IIFraction[sphere]);
    /* 2018.02.01 added */

    /* 2018.03.18 added */
    if (sscanf(line, "CollapseTestSphereNumberOfBubbles[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereNumberOfBubbles[%d] = %d", &sphere,
          &NumberOfBubbles[sphere]);
    if (sscanf(line, "CollapseTestSphereMetallicityDistributionCase[%d]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereMetallicityDistributionCase[%d] = %d", &sphere,
          &MetallicityDistributionCase[sphere]);  
    /* 2018.03.18 added */

    if (sscanf(line, "CollapseTestSphereVelocityTowardCenter[%d]", &sphere) > 0)
    ret += sscanf(line, "CollapseTestSphereVelocityTowardCenter[%d] = %d", &sphere,
      &CollapseTestSphereVelocityTowardCenter[sphere]);

    ret += sscanf(line, "CollapseTestSphereDarkMatterHaloMass = %d",
    &CollapseTestSphereDarkMatterHaloMass);

    ret += sscanf(line, "CollapseTestInitialFractionHII = %f",
      &CollapseTestInitialFractionHII);
    ret += sscanf(line, "CollapseTestInitialFractionHeII = %f",
      &CollapseTestInitialFractionHeII);
    ret += sscanf(line, "CollapseTestInitialFractionHeIII = %f",
      &CollapseTestInitialFractionHeIII);
    ret += sscanf(line, "CollapseTestInitialFractionHM = %f",
      &CollapseTestInitialFractionHM);
    ret += sscanf(line, "CollapseTestInitialFractionH2I = %f",
      &CollapseTestInitialFractionH2I);
    ret += sscanf(line, "CollapseTestInitialFractionH2II = %f",
      &CollapseTestInitialFractionH2II);

    ret += sscanf(line, "LengthUnits = %f",
      &LengthUnits);


    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CollapseTest") 
    && line[0] != '#')
    fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up grid */

  /********************************************************************************************************************/
  /* 2018.08.07 modified, 2018.03.23 added */

  float rr = 0., x = 0., y = 0., z = 0.;
  float xpos, ypos, zpos;

  // T_cmb ~ 60-100 K, so T_SNR ~ 600-1000 K

  /***  core collapse SNe parameters, based on 15 Msun star form ICG_2016 (Ken) ***/
  /* r1 < r2 */
  float ccsne_r1 = 6.171355e20; // 0.2 kpc in cm
  float ccsne_r2 = 1.54285e+21; // 0.5 kpc in cm
  float ccsne_d1 = 1.673e-27; // 1.e-3 1/cc
  float ccsne_d2 = 2. * 1.673e-26; // 2.e-2 1/cc
  float ccsne_Z1 = 1.e-2; // in Z_sun
  float ccsne_Z2 = 1.e-4; // in Z_sun
  float ccsne_T1 = 1.e3; // Kelvin
  float ccsne_T2 = 6.e2; // Kelvin
  /***  core collapse SNe parameters  ***/

  /***  hypernovae parameters, based on 60 Msun star from ICG_2016 (Ken)  ***/
  /* r1 < r2 */
  float hne_r1 = 1.54285e+21; // 0.5 kpc in cm
  float hne_r2 = 3.0857e+21; // 1 kpc in cm
  float hne_d1 = 2.5 * 1.673e-27; // 2.5.e-3 1/cc
  float hne_d2 = 2.5 * 1.673e-26; // 2.5e-2 1/cc
  float hne_Z1 = 1.e-2; // in Z_sun
  float hne_Z2 = 1.e-3; // in Z_sun
  float hne_T1 = 1.5e3; // Kelvin
  float hne_T2 = 8.e2; // Kelvin
  /***  hypernovae parameters  ***/

  /***  pair instability SNe parameters  ***/
  /* r1 < r2 */
  float pisne_r1 = 2.3142e21; // 0.75 kpc in cm
  float pisne_r2 = 4.6285e21; // 1.5 kpc in cm
  float pisne_d1 = 3. * 1.673e-27; // 3.e-3 1/cc
  float pisne_d2 = 3. * 1.673e-26; // 3.e-2 1/cc
  float pisne_Z1 = 2.e-2; // in Z_sun
  float pisne_Z2 = 2.e-3; // in Z_sun
  float pisne_T1 = 2.e3; // Kelvin
  float pisne_T2 = 1.e3; // Kelvin
  /***  pair instability SNe parameters  ***/


  /* initialization of spheres */
  for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++) {
    // printf("yoyoyoyo\n");
    printf("%d\n", CollapseTestSphereDarkMatterHaloMass);
    // for (int bubble = 0; bubble < NumberOfBubbles[sphere]; bubble ++) {
    //   printf("BubbleDensity: %f\n", BubbleDensity[sphere][bubble]);
    // }

    float outer_radius = (CollapseTestSphereSmoothSurface[sphere] == 1) ? 
          CollapseTestSphereSmoothRadius[sphere]*CollapseTestSphereRadius[sphere] : CollapseTestSphereRadius[sphere];
    // printf("outer_radius: %f\n", outer_radius);

    int count = 0;
    float left_edge;
    left_edge = CollapseTestSpherePosition[sphere][0] - outer_radius;
    // printf("left_edge %f\n", left_edge);

    /* Decide the center of the bubbles */
    if (CollapseTestSphereDarkMatterHaloMass == 9) { 
      printf("1e9\n");  
      while (count < NumberOfBubbles[sphere]) {
        x = ( float(rand())/float(RAND_MAX) ) * (outer_radius * 2.) + left_edge;
        y = ( float(rand())/float(RAND_MAX) ) * (outer_radius * 2.) + left_edge;
        z = ( float(rand())/float(RAND_MAX) ) * (outer_radius * 2.) + left_edge;
        
        if (x < left_edge || x > left_edge+outer_radius*2.) 
          printf("OutofBoundary x: %f\n", x);
        else if (y < left_edge || y > left_edge+outer_radius*2.)
          printf("OutofBoundary y: %f\n", y);
        else if (z < left_edge || z > left_edge+outer_radius*2.)
          printf("OutofBoundary z: %f\n", z);
        
        xpos = x-CollapseTestSpherePosition[sphere][0];
        ypos = y-CollapseTestSpherePosition[sphere][1];
        zpos = z-CollapseTestSpherePosition[sphere][2];
        rr = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);

        float tx = x, ty = y, tz = z;
        if ( MetallicityDistributionCase[sphere] < 5) {
          if (rr <= 0.9 * outer_radius && count < NumberOfBubbles[sphere]-1) {
            center_of_bubble[sphere][count][0] = tx;
            center_of_bubble[sphere][count][1] = ty;
            center_of_bubble[sphere][count][2] = tz;
            // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
            count ++;
          }

          if (count == NumberOfBubbles[sphere]-1) {
            if ( abs(rr - 1.5428e21/LengthUnits) / (1.5428e21/LengthUnits) < 0.0001 ) {
              printf("rr: %.6e\n", rr);
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              count ++;
            }
          }
        }
        if ( MetallicityDistributionCase[sphere] == 5) {
          if ( (rr > 0.3) && (rr <= 0.9 * outer_radius) ) {
            center_of_bubble[sphere][count][0] = tx;
            center_of_bubble[sphere][count][1] = ty;
            center_of_bubble[sphere][count][2] = tz;
            // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
            count ++;
          }
        }
        if (MetallicityDistributionCase[sphere] == 6) {
          if (rr <= 0.6 * outer_radius) { // && count < NumberOfBubbles[sphere]-1
            center_of_bubble[sphere][count][0] = tx;
            center_of_bubble[sphere][count][1] = ty;
            center_of_bubble[sphere][count][2] = tz;
            // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
            count ++;
          }

          /*
          if (count == NumberOfBubbles[sphere]-1) {
              center_of_bubble[sphere][count][0] = 0.5;
              center_of_bubble[sphere][count][1] = 0.5;
              center_of_bubble[sphere][count][2] = 0.5;
              count ++;
          }
          */
        }

        if (MetallicityDistributionCase[sphere] == 7) {
          if (count < 5) {
            if ( rr <= 0.75 * outer_radius && ty > 0.5) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            }
          } else {
             if ( rr <= 0.75 * outer_radius && ty < 0.5) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            } 
          }
        } // end Shingo
            
        if (MetallicityDistributionCase[sphere] == 8) {
          if (count < 8) {
            if ( rr <= 0.9 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            }
          } else if (count > 8 && count < 16) {
            if ( rr <= 0.7 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            } 
          } else {
            if ( rr <= 0.5 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            }
          }
        } // end Salpeter

        if (MetallicityDistributionCase[sphere] == 9 ) {
          if ( rr <= 0.8 * outer_radius ) {
            center_of_bubble[sphere][count][0] = tx;
            center_of_bubble[sphere][count][1] = ty;
            center_of_bubble[sphere][count][2] = tz;
            // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
            count ++;
          }
        } // end Flat IMF
      } // end while
    } // 1e9 M_sun halo 
    else if (CollapseTestSphereDarkMatterHaloMass == 8) {
      printf("1e8\n");
      while (count < NumberOfBubbles[sphere]) {
        x = ( float(rand())/float(RAND_MAX) ) * (outer_radius * 2.) + left_edge;
        y = ( float(rand())/float(RAND_MAX) ) * (outer_radius * 2.) + left_edge;
        z = ( float(rand())/float(RAND_MAX) ) * (outer_radius * 2.) + left_edge;
        
        if (x < left_edge || x > left_edge+outer_radius*2.) 
          printf("OutofBoundary x: %f\n", x);
        else if (y < left_edge || y > left_edge+outer_radius*2.)
          printf("OutofBoundary y: %f\n", y);
        else if (z < left_edge || z > left_edge+outer_radius*2.)
          printf("OutofBoundary z: %f\n", z);
        
        xpos = x-CollapseTestSpherePosition[sphere][0];
        ypos = y-CollapseTestSpherePosition[sphere][1];
        zpos = z-CollapseTestSpherePosition[sphere][2];
        rr = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);

        float tx = x, ty = y, tz = z;
        if (MetallicityDistributionCase[sphere] == 7) { // Shingo, 3 SNR
          if (count < 3) {
            if (rr <= 0.5 * outer_radius && rr >= 0.3 * outer_radius) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f, rr: %f\n\n", count, x, y, z, rr);
              count ++;
            }
          }
        } // Shingo         
            
        if (MetallicityDistributionCase[sphere] == 8) { // Salpeter, 5 SNR
          if (count < 4) {
            if ( rr <= 0.9 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            }
          } else if (count == 4 ) {
            if ( rr <= 0.7 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            } 
          } 
        } // Salpeter

        if (MetallicityDistributionCase[sphere] == 9 ) { // Flat IMF, 3 SNR 
          if (count == 0) {
            if ( rr <= 0.85 * outer_radius && rr > 0.75 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            }
          } else if (count == 1) {
            if ( rr <= 0.75 * outer_radius && rr >= 0.6 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            }
          } else { 
            if ( rr <= 0.6 * outer_radius ) {
              center_of_bubble[sphere][count][0] = tx;
              center_of_bubble[sphere][count][1] = ty;
              center_of_bubble[sphere][count][2] = tz;
              // printf("Bubble[%d]Center           x: %f,  y: %f,  z: %f\n", count, x, y, z);
              count ++;
            }
          }  
        } // Flat IMF
      } // end while
    } // 1e8 M_sun halo
    /* END Decide the center of the bubbles */

    count = 0;
    /* Profile of the density/metallicity/temperature of each bubbles for case <= 6*/
    while (count < NumberOfBubbles[sphere]) {
      const float vx = 0.01 * float(rand())/float(RAND_MAX); 
      const float vy = 0.01 * float(rand())/float(RAND_MAX);
      const float vz = 0.01 * float(rand())/float(RAND_MAX);
      BubbleVelocity[sphere][count][0] = 0.; // different from the drift velocity toward the center, which is calculated in 
      BubbleVelocity[sphere][count][1] = 0.; // Grid_CollapseTestInitializeGrid.C
      BubbleVelocity[sphere][count][2] = 0.; 

      if (MetallicityDistributionCase[sphere] == 1 || MetallicityDistributionCase[sphere] == 4 || MetallicityDistributionCase[sphere] == 5) {
        BubbleRadius[sphere][count] = 3.0857e+21; // 1 kpc in cm
        BubbleInnerRadius[sphere][count] = 6.171355e20; // 200 pc in cm
        BubbleInnerDensity[sphere][count] = 1.673e-27; // 1.e-3 1/cc
        BubbleOuterDensity[sphere][count] = 2. * 1.673e-26; // 2.e-2 1/cc
        BubbleInnerMetallicity[sphere][count] = 1.e-2; // in Z_sun
        BubbleOuterMetallicity[sphere][count] = 1.e-4; // in Z_sun
        BubbleInnerTemperature[sphere][count] = 1.e4; // Kelvin
        BubbleOuterTemperature[sphere][count] = 5.e3; // Kelvin
        // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", count, vx, vy, vz);
        count ++;
      }
      if (MetallicityDistributionCase[sphere] == 2) {
        BubbleVelocity[sphere][count][0] = vx;
        BubbleVelocity[sphere][count][1] = vy;
        BubbleVelocity[sphere][count][2] = vz;

        BubbleDensity[sphere][count] = 1.e4;
        BubbleMetallicity[sphere][count] = 1.e-4;
        BubbleTemperature[sphere][count] = 1.e4;
        // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", count, vx, vy, vz);
        count ++;
      } 
      if (MetallicityDistributionCase[sphere] == 3) {
        BubbleVelocity[sphere][count][0] = vx;
        BubbleVelocity[sphere][count][1] = vy;
        BubbleVelocity[sphere][count][2] = vz;

        BubbleDensity[sphere][count] = 1.e4;
        BubbleMetallicity[sphere][count] = 1.e-4;
        BubbleTemperature[sphere][count] = 1.e4;
        count ++;
      }
      if (MetallicityDistributionCase[sphere] == 6) {
        BubbleRadius[sphere][count] = 3.0857e+21; // 1 kpc in cm
        BubbleInnerRadius[sphere][count] = 6.171355e20; // 200 pc in cm
        BubbleInnerDensity[sphere][count] = 1.673e-27; // 1.e-3 1/cc
        BubbleOuterDensity[sphere][count] = 2. * 1.673e-26; // 2.e-2 1/cc
        BubbleInnerMetallicity[sphere][count] = 1.e-2; // in Z_sun
        BubbleOuterMetallicity[sphere][count] = 1.e-4; // in Z_sun
        BubbleInnerTemperature[sphere][count] = 1.e4; // Kelvin
        BubbleOuterTemperature[sphere][count] = 5.e3; // Kelvin
        // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", count, vx, vy, vz);
        count ++;
      }
      if (MetallicityDistributionCase[sphere] > 6) {
        count = NumberOfBubbles[sphere];
        printf("HEYYOOO\n");
      }
      if (MetallicityDistributionCase[sphere] == -1){
        printf("No metal bubbles inside sphere: %d\n", sphere);
        break;
      }
    } /* END Profile of the density/metallicity/temperature of each bubbles for case <= 6*/


    /* Profile of the density/metallicity/temperature of each bubbles for case >= 7, 1e9 Msun Halo */
    if (CollapseTestSphereDarkMatterHaloMass == 9) {
      printf("1e9\n");
      for (int bubble=0; bubble<NumberOfBubbles[sphere]; bubble++) {
        // printf("aldkfjd\n");
        if (MetallicityDistributionCase[sphere] == 7) { // Shingo's first stars IMF, 10 SNR
          if (bubble < 2) {
            BubbleRadius[sphere][bubble] = ccsne_r2;
            BubbleInnerRadius[sphere][bubble] = ccsne_r1;
            BubbleInnerDensity[sphere][bubble] = ccsne_d1;
            BubbleOuterDensity[sphere][bubble] = ccsne_d2;
            BubbleInnerMetallicity[sphere][bubble] = ccsne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = ccsne_Z2;
            BubbleInnerTemperature[sphere][bubble] = ccsne_T1;
            BubbleOuterTemperature[sphere][bubble] = ccsne_T2;
          } else {
            BubbleRadius[sphere][bubble] = pisne_r2;
            BubbleInnerRadius[sphere][bubble] = pisne_r1;
            BubbleInnerDensity[sphere][bubble] = pisne_d1;
            BubbleOuterDensity[sphere][bubble] = pisne_d2;
            BubbleInnerMetallicity[sphere][bubble] = pisne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = pisne_Z2;
            BubbleInnerTemperature[sphere][bubble] = pisne_T1;
            BubbleOuterTemperature[sphere][bubble] = pisne_T2;
          }
        } /* END case 7 */

        if (MetallicityDistributionCase[sphere] == 8 ) { // Salpeter's IMF with peak at 10 Msun, 19 SNR
          if (bubble < 18) {
            BubbleRadius[sphere][bubble] = ccsne_r2;
            BubbleInnerRadius[sphere][bubble] = ccsne_r1;
            BubbleInnerDensity[sphere][bubble] = ccsne_d1;
            BubbleOuterDensity[sphere][bubble] = ccsne_d2;
            BubbleInnerMetallicity[sphere][bubble] = ccsne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = ccsne_Z2;
            BubbleInnerTemperature[sphere][bubble] = ccsne_T1;
            BubbleOuterTemperature[sphere][bubble] = ccsne_T2;
          } else {
            BubbleRadius[sphere][bubble] = hne_r2;
            BubbleInnerRadius[sphere][bubble] = hne_r1;
            BubbleInnerDensity[sphere][bubble] = hne_d1;
            BubbleOuterDensity[sphere][bubble] = hne_d2;
            BubbleInnerMetallicity[sphere][bubble] = hne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = hne_Z2;
            BubbleInnerTemperature[sphere][bubble] = hne_T1;
            BubbleOuterTemperature[sphere][bubble] = hne_T2;
          }        
        } /* END case 8 */

        if (MetallicityDistributionCase[sphere] == 9) { // test, 3 ccsne, 3 hne, 3 pisne
          if (bubble < 3) {
            BubbleRadius[sphere][bubble] = ccsne_r2;
            BubbleInnerRadius[sphere][bubble] = ccsne_r1;
            BubbleInnerDensity[sphere][bubble] = ccsne_d1;
            BubbleOuterDensity[sphere][bubble] = ccsne_d2;
            BubbleInnerMetallicity[sphere][bubble] = ccsne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = ccsne_Z2;
            BubbleInnerTemperature[sphere][bubble] = ccsne_T1;
            BubbleOuterTemperature[sphere][bubble] = ccsne_T2;
            // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", bubble, vx, vy, vz);
          } else if (bubble >=3 && bubble < 6) {
            BubbleRadius[sphere][bubble] = hne_r2;
            BubbleInnerRadius[sphere][bubble] = hne_r1;
            BubbleInnerDensity[sphere][bubble] = hne_d1;
            BubbleOuterDensity[sphere][bubble] = hne_d2;
            BubbleInnerMetallicity[sphere][bubble] = hne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = hne_Z2;
            BubbleInnerTemperature[sphere][bubble] = hne_T1;
            BubbleOuterTemperature[sphere][bubble] = hne_T2;
            // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", bubble, vx, vy, vz);         
          } else if (bubble >= 6 && bubble < 9) {
            BubbleRadius[sphere][bubble] = pisne_r2;
            BubbleInnerRadius[sphere][bubble] = pisne_r1;
            BubbleInnerDensity[sphere][bubble] = pisne_d1;
            BubbleOuterDensity[sphere][bubble] = pisne_d2;
            BubbleInnerMetallicity[sphere][bubble] = pisne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = pisne_Z2;
            BubbleInnerTemperature[sphere][bubble] = pisne_T1;
            BubbleOuterTemperature[sphere][bubble] = pisne_T2;
            // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", bubble, vx, vy, vz);
          }
        } /* END case 9 */
      } /* END Profile of the density/metallicity/temperature of each bubbles for case >= 7, 1e9 Msun Halo */
    }


    /* Profile of the density/metallicity/temperature of each bubbles for case >= 7, 1e8 Msun Halo */
    if (CollapseTestSphereDarkMatterHaloMass == 8) {
      printf("1e8\n");
      printf("%d\n", NumberOfBubbles[sphere]);
      for (int bubble=0; bubble<NumberOfBubbles[sphere]; bubble++) {
        // printf("aldkfjd\n");
        if (MetallicityDistributionCase[sphere] == 7) { // Shingo's first stars IMF, 3 psne
          printf("Shingo\n");
          if (bubble < 0) {
            printf("bubble0\n");
            BubbleRadius[sphere][bubble] = ccsne_r2;
            BubbleInnerRadius[sphere][bubble] = ccsne_r1;
            BubbleInnerDensity[sphere][bubble] = ccsne_d1;
            BubbleOuterDensity[sphere][bubble] = ccsne_d2;
            BubbleInnerMetallicity[sphere][bubble] = ccsne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = ccsne_Z2;
            BubbleInnerTemperature[sphere][bubble] = ccsne_T1;
            BubbleOuterTemperature[sphere][bubble] = ccsne_T2;
          } else {
            printf("testttt\n");
            BubbleRadius[sphere][bubble] = pisne_r2;
            BubbleInnerRadius[sphere][bubble] = pisne_r1;
            BubbleInnerDensity[sphere][bubble] = pisne_d1;
            BubbleOuterDensity[sphere][bubble] = pisne_d2;
            BubbleInnerMetallicity[sphere][bubble] = pisne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = pisne_Z2;
            BubbleInnerTemperature[sphere][bubble] = pisne_T1;
            BubbleOuterTemperature[sphere][bubble] = pisne_T2;
          }
        } /* END case 7 */

        if (MetallicityDistributionCase[sphere] == 8 ) { // Salpeter's IMF with peak at 10 Msun, 5 SNR
          if (bubble < 4) {
            BubbleRadius[sphere][bubble] = ccsne_r2;
            BubbleInnerRadius[sphere][bubble] = ccsne_r1;
            BubbleInnerDensity[sphere][bubble] = ccsne_d1;
            BubbleOuterDensity[sphere][bubble] = ccsne_d2;
            BubbleInnerMetallicity[sphere][bubble] = ccsne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = ccsne_Z2;
            BubbleInnerTemperature[sphere][bubble] = ccsne_T1;
            BubbleOuterTemperature[sphere][bubble] = ccsne_T2;
          } else {
            BubbleRadius[sphere][bubble] = hne_r2;
            BubbleInnerRadius[sphere][bubble] = hne_r1;
            BubbleInnerDensity[sphere][bubble] = hne_d1;
            BubbleOuterDensity[sphere][bubble] = hne_d2;
            BubbleInnerMetallicity[sphere][bubble] = hne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = hne_Z2;
            BubbleInnerTemperature[sphere][bubble] = hne_T1;
            BubbleOuterTemperature[sphere][bubble] = hne_T2;
          }        
        } /* END case 8 */

        if (MetallicityDistributionCase[sphere] == 9) { // Flat IMF, 1 ccsne, 1 hne, 1 pisne
          if (bubble == 0) {
            BubbleRadius[sphere][bubble] = ccsne_r2;
            BubbleInnerRadius[sphere][bubble] = ccsne_r1;
            BubbleInnerDensity[sphere][bubble] = ccsne_d1;
            BubbleOuterDensity[sphere][bubble] = ccsne_d2;
            BubbleInnerMetallicity[sphere][bubble] = ccsne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = ccsne_Z2;
            BubbleInnerTemperature[sphere][bubble] = ccsne_T1;
            BubbleOuterTemperature[sphere][bubble] = ccsne_T2;
            // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", bubble, vx, vy, vz);
          } else if (bubble == 1) {
            BubbleRadius[sphere][bubble] = hne_r2;
            BubbleInnerRadius[sphere][bubble] = hne_r1;
            BubbleInnerDensity[sphere][bubble] = hne_d1;
            BubbleOuterDensity[sphere][bubble] = hne_d2;
            BubbleInnerMetallicity[sphere][bubble] = hne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = hne_Z2;
            BubbleInnerTemperature[sphere][bubble] = hne_T1;
            BubbleOuterTemperature[sphere][bubble] = hne_T2;
            // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", bubble, vx, vy, vz);         
          } else if (bubble == 2) {
            BubbleRadius[sphere][bubble] = pisne_r2;
            BubbleInnerRadius[sphere][bubble] = pisne_r1;
            BubbleInnerDensity[sphere][bubble] = pisne_d1;
            BubbleOuterDensity[sphere][bubble] = pisne_d2;
            BubbleInnerMetallicity[sphere][bubble] = pisne_Z1;
            BubbleOuterMetallicity[sphere][bubble] = pisne_Z2;
            BubbleInnerTemperature[sphere][bubble] = pisne_T1;
            BubbleOuterTemperature[sphere][bubble] = pisne_T2;
            // printf("Bubble[%d]CenterVelocity  vx: %f, vy: %f, vz: %f\n", bubble, vx, vy, vz);
          }
        } /* END case 9 */
      } /* END Profile of the density/metallicity/temperature of each bubbles for case >= 7, 1e9 Msun Halo */
    }


  } /* END initialization of spheres */

  printf("before output\n");

  FILE *Outfptr;
  char file_name_[MAX_LINE_LENGTH];
  

  for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++) {
    /* create txt file to write */
    for (int mkk = 1; mkk < 10; mkk++) {
      if (MetallicityDistributionCase[sphere] == mkk) {
        sprintf(file_name_, "CTNC_BubbleInitialize_sphere%d_metalic%d.txt", sphere, mkk);
        printf("%s\n", file_name_);
      } 
    }

    Outfptr = fopen(file_name_, "w");
    if (MetallicityDistributionCase[sphere] == 1 || MetallicityDistributionCase[sphere] == 4 || MetallicityDistributionCase[sphere] == 5) {
       for (bubble = 0; bubble < NumberOfBubbles[sphere]; bubble++) {
        fprintf(Outfptr, "CollapseTestSphereBubbleCenter[%d][%d] = %f %f %f\n",
                sphere, bubble, center_of_bubble[sphere][bubble][0],
                center_of_bubble[sphere][bubble][1],
                center_of_bubble[sphere][bubble][2]);

        fprintf(Outfptr, "CollapseTestSphereBubbleInnerDensity[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerDensity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterDensity[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterDensity[sphere][bubble]);
        
        fprintf(Outfptr, "CollapseTestSphereBubbleInnerMetallicity[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerMetallicity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterMetallicity[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterMetallicity[sphere][bubble]);

        fprintf(Outfptr, "CollapseTestSphereBubbleInnerTemperature[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerTemperature[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterTemperature[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterTemperature[sphere][bubble]);
        
        fprintf(Outfptr, "CollapseTestSphereBubbleRadius[%d][%d] = %g\n",
                sphere, bubble, BubbleRadius[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleInnerRadius[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerRadius[sphere][bubble]);
      }     
    } else if ( MetallicityDistributionCase[sphere] == 6) {
       for (bubble = 0; bubble < NumberOfBubbles[sphere]; bubble++) {
        fprintf(Outfptr, "CollapseTestSphereBubbleCenter[%d][%d] = %f %f %f\n",
                sphere, bubble, center_of_bubble[sphere][bubble][0],
                center_of_bubble[sphere][bubble][1],
                center_of_bubble[sphere][bubble][2]);

        fprintf(Outfptr, "CollapseTestSphereBubbleInnerDensity[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerDensity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterDensity[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterDensity[sphere][bubble]);
        
        fprintf(Outfptr, "CollapseTestSphereBubbleInnerMetallicity[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerMetallicity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterMetallicity[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterMetallicity[sphere][bubble]);

        fprintf(Outfptr, "CollapseTestSphereBubbleInnerTemperature[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerTemperature[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterTemperature[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterTemperature[sphere][bubble]);
        
        fprintf(Outfptr, "CollapseTestSphereBubbleRadius[%d][%d] = %g\n",
                sphere, bubble, BubbleRadius[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleInnerRadius[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerRadius[sphere][bubble]);
      }     
    } else if (MetallicityDistributionCase[sphere] == 2 || MetallicityDistributionCase[sphere] == 3) {
      for (bubble = 0; bubble < NumberOfBubbles[sphere]; bubble++) {
        fprintf(Outfptr, "CollapseTestSphereBubbleCenter[%d][%d] = %f %f %f\n",
                sphere, bubble, center_of_bubble[sphere][bubble][0],
                center_of_bubble[sphere][bubble][1],
                center_of_bubble[sphere][bubble][2]);
        fprintf(Outfptr, "CollapseTestSphereBubbleVelocity[%d][%d] = %f %f %f\n",
                sphere, bubble, BubbleVelocity[sphere][bubble][0],
                BubbleVelocity[sphere][bubble][1],
                BubbleVelocity[sphere][bubble][2]);
        fprintf(Outfptr, "CollapseTestSphereBubbleDensity[%d][%d] = %g\n",
                sphere, bubble, BubbleDensity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleMetallicity[%d][%d] = %g\n",
                sphere, bubble, BubbleMetallicity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleTemperature[%d][%d] = %g\n",
                sphere, bubble, BubbleTemperature[sphere][bubble]);
      }      
    } else if ( MetallicityDistributionCase[sphere] == 7 || MetallicityDistributionCase[sphere] == 8 || MetallicityDistributionCase[sphere] == 9 ) {
       for (bubble = 0; bubble < NumberOfBubbles[sphere]; bubble++) {
        fprintf(Outfptr, "CollapseTestSphereBubbleCenter[%d][%d] = %f %f %f\n",
                sphere, bubble, center_of_bubble[sphere][bubble][0],
                center_of_bubble[sphere][bubble][1],
                center_of_bubble[sphere][bubble][2]);

        fprintf(Outfptr, "CollapseTestSphereBubbleInnerDensity[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerDensity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterDensity[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterDensity[sphere][bubble]);

        fprintf(Outfptr, "CollapseTestSphereBubbleInnerMetallicity[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerMetallicity[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterMetallicity[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterMetallicity[sphere][bubble]);

        fprintf(Outfptr, "CollapseTestSphereBubbleInnerTemperature[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerTemperature[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleOuterTemperature[%d][%d] = %g\n",
                sphere, bubble, BubbleOuterTemperature[sphere][bubble]);
        
        fprintf(Outfptr, "CollapseTestSphereBubbleRadius[%d][%d] = %g\n",
                sphere, bubble, BubbleRadius[sphere][bubble]);
        fprintf(Outfptr, "CollapseTestSphereBubbleInnerRadius[%d][%d] = %g\n",
                sphere, bubble, BubbleInnerRadius[sphere][bubble]);
      }     
    }

    fclose(Outfptr);

  }

  printf("Writing Bubble Info\n");

  return 0;
  /* 2018.08.07 modified, 2018.03.23 added */
  /********************************************************************************************************************/

}
