/***********************************************************************
/
/  POP III IMF LOOKUP TABLE INITIALIZATION
/
/  written by: John Wise
/  date:       April, 2010
/  modified1:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

void mt_init(unsigned_int seed);
unsigned_long_int mt_random(void);

int StarParticlePopIII_ShingoIMFInitialize(void)
{
  //fprintf(stdout, "Shingo PopIIIInitialMassFunction\n");
  const float CutoffExponent = 1.6;

  if (IMFData != NULL)
    return SUCCESS;

  IMFData = new float[IMF_TABLE_ENTRIES];

  int i;
  float m, m0, dm, total_fn;

  dm = log10(PopIIIUpperMassCutoff / PopIIILowerMassCutoff) / 
    (float) (IMF_TABLE_ENTRIES-1);
  m0 = log10(PopIIILowerMassCutoff);
  total_fn = 0;

  float slope_1 = -1.*(2-0.477)/(2.813-2.398);
  float slope_2 = -1.*(2-0.602)/(2.398-1.845);
  float slope_3 = -1.*(log10(40)-log10(4))/(log10(70)-log10(25));
  float slope_4 = -1.*(log10(40)-log10(6))/(log10(25)-log10(10));
  float MStar_1 = 250;
  float MStar_2 = 70;
  float MStar_3 = 25;

  float coeff_a = POW(MStar_1, -2.528) * POW(MStar_3, -2.236);
  float coeff_b = POW(MStar_2, -2.528-2.236);

  for (i = 0; i < IMF_TABLE_ENTRIES; i++) {
      m = POW(10.0, m0 + i*dm);
      if (PopIIIStarMass < m && m <= PopIIIUpperMassCutoff)
        total_fn += POW(m/MStar_1, -3.67);
      else if (70 < m && m <= PopIIIStarMass)
        total_fn += POW(MStar_1/m, -2.528);
      else if (25 < m && m <= 70)
        total_fn += (coeff_a/coeff_b) * POW(m/MStar_3, -2.236);
      else if (10 < m && m <= 25)
        total_fn += (coeff_a/coeff_b) * POW(MStar_3/m, -2.07);
      IMFData[i] = total_fn;
  }

  // Normalize
  for (i = 0; i < IMF_TABLE_ENTRIES; i++)
    IMFData[i] /= IMFData[IMF_TABLE_ENTRIES-1];

  /* Initialize the random number generator.  Call it 1+NumberOfCalls
     to get the generator back to the same state as it was before the
     restart (if any). */

  if (PopIIIInitialMassFunctionSeed == INT_UNDEFINED)
    mt_init(time(NULL)); //+100*MyProcessorNumber);
  else
    mt_init(PopIIIInitialMassFunctionSeed); //+100*MyProcessorNumber);
  unsigned_long_int trash;
  for (i = 0; i < 1+PopIIIInitialMassFunctionCalls; i++)
    trash = mt_random();

  return SUCCESS;

}
