	  i%  C   k820309              15.0        &W                                                                                                           
       ../src/fourier_taylor_tools.f90 FOURIER_TAYLOR_TOOLS                                                    
                                                        
                 
          +       -DTû!@        6.283185307179586476925286766559005768394D0#         @                                                      #CALCULATE_FOURIER_TAYLOR_BIASES%TRIM    #CALCULATE_FOURIER_TAYLOR_BIASES%SQRT    #CALCULATE_FOURIER_TAYLOR_BIASES%INT    #CALCULATE_FOURIER_TAYLOR_BIASES%DBLE    #N    #NORDER 	   #VERBOSE 
   #FILEOUT    #RESIDUP1    #FACTO    #ETA                                                   TRIM                                                SQRT                                                INT                                                DBLE           D @                                                     D @                               	                                                       
                       @                                                   1          D                                                    
     p           & p          5  p        r    n                                       2       5  p        r    n                                      2p         p                                            D @                                                  
     p           & p         5  p        r 	         5  p        r 	   p         p                                   D @                                                  
       p            5  p        r    n                                           2  5  p        r    n                                          2p          5  p        r    n                                      2  &   5  p        r    n                                          2 5  p        r    n                                      2  & p         5  p        r 	          5  p        r    n                                      2  5  p        r    n                                          2p            5  p        r 	   p         p                                            #         @                                                     #COMPUTE_FACTO%DBLE    #NORDER    #FACTO                                                   DBLE                                                                 D                                                    
     p           & p         5  p        r          5  p        r    p         p                          #         @                                                     #COMPUTE_ETA%COS    #COMPUTE_ETA%SIN    #COMPUTE_ETA%DBLE    #KNY    #NORDER    #ETA                                                   COS                                                SIN                                                DBLE            @                                                                                                           D                                                    
       p           5  p        r     5  p        r    p         5  p        r      &  5  p        r    5  p        r      & p         5  p        r          5  p        r     5  p        r    p            5  p        r    p         p                          %         @                                                  
       #SHOTNOISE_FACTOR%MOD    #I    #J    #K    #NORDER    #N     #FACTO !                                                  MOD           D @                                                     D @                                                     D @                                                                                                            D @                                                                                     !                    
 
    p           & p         5  p        r          5  p        r    p         p                          %         @                               "                   
       #SHOTNOISE_FACTOR_ASYMP%MOD #   #SHOTNOISE_FACTOR_ASYMP%DBLE $   #I %   #J &   #K '   #NORDER (   #N )   #FACTO *                                             #     MOD                                           $     DBLE           D @                               %                      D @                               &                      D @                               '                       @                               (                      D @                               )                                                     *                    
 	    p           & p         5  p        r (         5  p        r (   p         p                          %         @                               +                    
       #I ,   #J -   #K .   #NORDER /   #N 0   #FACTO 1   #ETA 2             D @                               ,                      D @                               -                      D @                               .                                                       /                      D @                               0                     D @                              1                    
     p           & p         5  p        r /         5  p        r /   p         p                                   D @                              2                    
       p            5  p        r 0   n                                           2  5  p        r 0   n                                          2p          5  p        r 0   n                                      2  &   5  p        r 0   n                                          2 5  p        r 0   n                                      2  & p         5  p        r /          5  p        r 0   n                                      2  5  p        r 0   n                                          2p            5  p        r /   p         p                                            %         @                               3                   
       #DELTA_MOM%DBLE 4   #I 5   #J 6   #K 7   #NORDER 8   #N 9                                             4     DBLE            @                               5                       @                               6                       @                               7                                                       8                       @                               9            #         @                                  :                    #I ;   #J <   #K =   #NORDER >   #KAPPA ?   #FACTO @   #ETA A   #N B                                              ;                                                       <                                                       =                                                       >                      D                                ?     
                                                @                    
     p           & p         5  p        r >         5  p        r >   p         p                                                                   A                    
       p            5  p        r B   n                                           2  5  p        r B   n                                          2p          5  p        r B   n                                      2  &   5  p        r B   n                                          2 5  p        r B   n                                      2  & p         5  p        r >          5  p        r B   n                                      2  5  p        r B   n                                          2p            5  p        r >   p         p                                                                                       B                   =      fn#fn    Ý   @   J   TWOPIDEF             TWOPI+TWOPIDEF 0   ¸  >      CALCULATE_FOURIER_TAYLOR_BIASES 5   ö  =      CALCULATE_FOURIER_TAYLOR_BIASES%TRIM 5   3  =      CALCULATE_FOURIER_TAYLOR_BIASES%SQRT 4   p  <      CALCULATE_FOURIER_TAYLOR_BIASES%INT 5   ¬  =      CALCULATE_FOURIER_TAYLOR_BIASES%DBLE 2   é  @   a   CALCULATE_FOURIER_TAYLOR_BIASES%N 7   )  @   a   CALCULATE_FOURIER_TAYLOR_BIASES%NORDER 8   i  @   a   CALCULATE_FOURIER_TAYLOR_BIASES%VERBOSE 8   ©  L   a   CALCULATE_FOURIER_TAYLOR_BIASES%FILEOUT 9   õ  V  a   CALCULATE_FOURIER_TAYLOR_BIASES%RESIDUP1 6   K  ä   a   CALCULATE_FOURIER_TAYLOR_BIASES%FACTO 4   /    a   CALCULATE_FOURIER_TAYLOR_BIASES%ETA    ²
  w       COMPUTE_FACTO #   )  =      COMPUTE_FACTO%DBLE %   f  @   a   COMPUTE_FACTO%NORDER $   ¦  ä   a   COMPUTE_FACTO%FACTO      ¦       COMPUTE_ETA     0  <      COMPUTE_ETA%COS     l  <      COMPUTE_ETA%SIN !   ¨  =      COMPUTE_ETA%DBLE     å  @   a   COMPUTE_ETA%KNY #   %  @   a   COMPUTE_ETA%NORDER     e  ô  a   COMPUTE_ETA%ETA !   Y         SHOTNOISE_FACTOR %   ö  <      SHOTNOISE_FACTOR%MOD #   2  @   a   SHOTNOISE_FACTOR%I #   r  @   a   SHOTNOISE_FACTOR%J #   ²  @   a   SHOTNOISE_FACTOR%K (   ò  @   a   SHOTNOISE_FACTOR%NORDER #   2  @   a   SHOTNOISE_FACTOR%N '   r  ä   a   SHOTNOISE_FACTOR%FACTO '   V  Ä       SHOTNOISE_FACTOR_ASYMP +     <      SHOTNOISE_FACTOR_ASYMP%MOD ,   V  =      SHOTNOISE_FACTOR_ASYMP%DBLE )     @   a   SHOTNOISE_FACTOR_ASYMP%I )   Ó  @   a   SHOTNOISE_FACTOR_ASYMP%J )     @   a   SHOTNOISE_FACTOR_ASYMP%K .   S  @   a   SHOTNOISE_FACTOR_ASYMP%NORDER )     @   a   SHOTNOISE_FACTOR_ASYMP%N -   Ó  ä   a   SHOTNOISE_FACTOR_ASYMP%FACTO    ·         POWER_FACTOR    C  @   a   POWER_FACTOR%I      @   a   POWER_FACTOR%J    Ã  @   a   POWER_FACTOR%K $     @   a   POWER_FACTOR%NORDER    C  @   a   POWER_FACTOR%N #     ä   a   POWER_FACTOR%FACTO !   g    a   POWER_FACTOR%ETA    ê         DELTA_MOM    v  =      DELTA_MOM%DBLE    ³  @   a   DELTA_MOM%I    ó  @   a   DELTA_MOM%J    3  @   a   DELTA_MOM%K !   s  @   a   DELTA_MOM%NORDER    ³  @   a   DELTA_MOM%N    ó         COMPUTE_KAPPA       @   a   COMPUTE_KAPPA%I     Â  @   a   COMPUTE_KAPPA%J        @   a   COMPUTE_KAPPA%K %   B   @   a   COMPUTE_KAPPA%NORDER $      @   a   COMPUTE_KAPPA%KAPPA $   Â   ä   a   COMPUTE_KAPPA%FACTO "   ¦!    a   COMPUTE_KAPPA%ETA     )%  @   a   COMPUTE_KAPPA%N 