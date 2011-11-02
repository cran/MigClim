/*
** barriers.c: Functions for performing the barriers related steps.
**
** Wim Hordijk    Last modified: 03 October 2011
**
** This C code is based on the original Visual Basic code of Robin Engler.
*/

#include "migclim.h"


/*
** mcFilterByBarrier: Filter the current state matrix using the barrier matrix.
**                    Every occupied pixel that is part of a barrier is reset
**                    to state 0.
**
** Parameters:
**   - curState: A pointer to the current state matrix.
**   - barriers: A pointer to the barriers matrix.
*/

void mcFilterByBarrier (int **curState, int **barriers)
{
  int i, j;

  /*
  ** Filter the current state matrix.
  */
  for (i = 0; i < nrRows; i++)
  {
    for (j = 0; j < nrCols; j++)
    {
      if (barriers[i][j] > 0)
      {
	curState[i][j] = 0;
      }
    }
  }
}


/*
** mcIntersectsBarrier: Check whether there is a barrier between the source
**                      and sink pixels.
** Parameters:
**   - snkX:     The x-coordinate of the sink pixel.
**   - snkY:     The y-coordinate of the sink pixel.
**   - srcX:     The x-coordinate of the source pixel.
**   - srcY:     The y-coordinate of the source pixel.
**   - barriers: A pointer to the barriers matrix.
**
** Returns:
**   If there is a barrier: True.
**   Otherwise:             False.
*/

bool mcIntersectsBarrier (int snkX, int snkY, int srcX, int srcY,
			  int **barriers)
{
  int  dstX, dstY, i, pxlX, pxlY, distMax, barCounter;
  bool barFound;

  barFound = false;
  
  /*
  ** Calculate the distance in both dimensions between the source and sink
  ** pixels and take the largest of the two.
  */
  dstX = srcX - snkX;
  dstY = srcY - snkY;
  if (abs (dstX) >= abs (dstY))
  {
    distMax = abs (dstX);
  }
  else
  {
    distMax = abs (dstY);
  }

  /*
  ** Check the possible paths from source to sink and see if there is a path
  ** without barriers.
  */
  if (barrierType == WEAK_BARRIER)
  {
    /*
    ** Weak barrier: If there is at least one free path we're good.
    **
    ** BARRIER MIDDLE
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] > 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** BARRIER TOP_LEFT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY - 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] > 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** BARRIER TOP_RIGHT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY - 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] > 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** Barrier DOWN_LEFT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] > 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
    /*
    ** Barrier DOWN_RIGHT
    */
    barFound = false;
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + 0.49 + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] > 1)
      {
	barFound = true;
	break;
      }
    }
    if (!barFound)
    {
      goto End_of_Routine;
    }
  }
  else if (barrierType == STRONG_BARRIER)
  {
    /*
    ** Strong barrier: If more than one way is blocked by a barrier then
    **                 colonization fails.
    */
    barCounter = 0;
    barFound = false;
    /*
    ** BARRIER MIDDLE
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + (1.0 * i / distMax * dstX));
      pxlY = (int)round (snkY + (1.0 * i / distMax * dstY));
      if (barriers[pxlX][pxlY] > 1)
      {
	barCounter++;
	break;
      }
    }
    /*
    ** BARRIER TOP_LEFT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY - 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] > 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
    /*    
    ** BARRIER TOP_RIGHT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY - 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] > 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
    /*
    ** BARRIER DOWN_LEFT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX - 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY + 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] > 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
    /*
    ** BARRIER DOWN_RIGHT
    */
    for (i = 1; i <= distMax; i++)
    {
      pxlX = (int)round (snkX + 0.49 + (((i-1.0) / distMax * dstX) +
					((1.0 / distMax * dstX) / 2.0)));
      pxlY = (int)round (snkY + 0.49 + (((i-1.0) / distMax * dstY) +
					((1.0 / distMax * dstY) / 2.0)));
      if (barriers[pxlX][pxlY] > 1)
      {
	barCounter++;
	break;
      }
    }
    if (barCounter > 1)
    {
      barFound = true;
      goto End_of_Routine;
    }
  }        

 End_of_Routine:
  /*
  ** Return the result.
  */
  return (barFound);
}
    

/*
** EoF: barriers.c
*/
