C========+=========+=========+=========+=========+=========+=========+=$
C   Mersenne Twister  (mt19937.f)
C
C     generate uniformly distributed random numbers between 0 and 1 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine sgrnd(seed)
*
      implicit integer(a-z)
*
* Period parameters
      parameter(N     =  624)
*
      dimension mt(0:N-1)
*                     the array for the state vector
      common /block/mti,mt
      save   /block/
*
*      setting initial seeds to mt[N] using
*      the generator Line 25 of Table 1 in
*      [KNUTH 1981, The Art of Computer Programming
*         Vol. 2 (2nd Ed.), pp102]
*
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue
*
      return
      end

************************************************************************
      double precision function grnd()
*
      implicit integer(a-z)
*
* Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
*                                    constant vector a
      parameter(UMASK = -2147483648)
*                                    most significant w-r bits
      parameter(LMASK =  2147483647)
*                                    least significant r bits
* Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
*
      dimension mt(0:N-1)
*                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
*                     mti==N+1 means mt[N] is not initialized
*
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
*                        mag01(x) = x * MATA for x=0,1
*
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
*
      if(mti.ge.N) then
*                       generate N words at one time
        if(mti.eq.N+1) then
*                            if sgrnd() has not been called,
          call sgrnd(4357)
*                              a default initial seed is used
        endif
*
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
*
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
*
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif
*
      return
      end

