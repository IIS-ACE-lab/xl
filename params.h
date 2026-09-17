#ifndef _DEFINE_H_
#define _DEFINE_H_

#ifdef OPEN_MPI
  #if QQ == 2
    static const unsigned BW_M = 512 * MPI_BLOCKS;
  #else // #if QQ == 2
    static const unsigned BW_M = 256;
  #endif // #if QQ == 2
#else // #ifdef OPEN_MPI
  #if QQ == 2
    static const unsigned BW_M = 512;
  #else // #if QQ == 2
    static const unsigned BW_M = 128;
  #endif // #if QQ == 2
#endif // #ifdef OPEN_MPI

   static const unsigned BW_N = 128;
   static const unsigned NSOL = 1;

#endif

