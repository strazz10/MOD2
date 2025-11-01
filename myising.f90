module parameters
  logical  :: periodic
  integer  :: lenght, dim_hilbert
  integer, dimension(:,:), allocatable :: spin_states
  real(8) :: g, lambda, lambdaAmul
end module parameters

program ising_quantum  		  !!test per modulo 2 metodi, ground state con davidson (almeno l'energia funziona)
     use parameters
     implicit none 
     integer :: i, j, itemp
     integer :: kmax          !!parametro davidson, autovalore massimo che la routine trova	
     real(8) :: cputime
     real(8), dimension(:), allocatable :: en, magn_long, magn_trasv      !!long=z, trasv=x
     double complex, dimension(:), allocatable :: ground_state   
     logical :: measure_magn, save_time  
    
     open(unit=10, file='parameters_in.dat', status='old')     !!leggo i parametri rilevanti
		read (10,*) lenght 		 !!lunghezza catena
		read (10,*) g 		     !!campo trasverso
		read (10,*) lambda 		 !!campo longitudinale
		read (10,*) periodic     !!condizioni al bordo periodiche se .true.
		read (10,*) measure_magn !!calcola o no magn.
	 close(10)
	 
	 dim_hilbert = 2**lenght				!!dimensione dello spazio di hilbert(= dim base computazionale)
	 allocate(spin_states(dim_hilbert, lenght))
	 spin_states = 0
	 
	 do i=1,dim_hilbert             !!inizializzo gli spin e li codifico in binario
		 itemp = i-1                !!per L =3, da 0 a 7 in binario, codifica al contrario del goat
		 do j=1,lenght
			 spin_states(i,lenght+1-j)=mod(itemp, 2)
			 itemp = itemp/2
		 end do
	 !print*, i,  spin_states(i, :)      !!verifica stati
	 end do
		 
	 allocate(ground_state(dim_hilbert), magn_long(lenght), magn_trasv(lenght))
	 kmax = 3
	 allocate(en(kmax))
	 lambdaAmul = lambda
	 
	 call davidson(dim_hilbert, kmax, en, ground_state)    !!magia: trovo il ground state (energia e funzione d'onda)
	 print *,'First three energy levels:   ', en(1:3)
	 
	 if (measure_magn) then								   !!calcolo di magnetizzazione per il ground state 
	     call total_magn(ground_state, magn_trasv, magn_long)   
	     if (periodic) then
	         print *,'GS magnetization along X:   ', magn_trasv(1)
	         print *,'GS magnetization along Z:   ', magn_long(1)
	     else
	         write (6,25)  magn_trasv
	         write (6,26)  magn_long     
	         25  format ('GS magnetization profile along X:   ', 100(es15.8,3x))
	         26  format ('GS magnetization profile along Z:   ', 100(es15.8,3x))
	     end if
	  end if 
	 
	 save_time = .true.  	!!printa il tempo o no
	 if (save_time) then
	     call cpu_time(cputime)
         print *,'Total CPU time:	', cputime
     end if
	 
	 stop               	!!dealloca e ferma tutto
  
end program ising_quantum 

subroutine total_magn(psi, magX, magZ)  
  use parameters  
  implicit none 
  real(8), dimension(lenght) :: magX, magZ
  real(8) :: sum_z
  double complex, dimension(dim_hilbert) :: psi
  double complex :: sum_x
  integer :: i, j, jtemp, temp, magn_temp
  
  sum_z = 0.d0               
  do i=1, dim_hilbert			!!magnetizzazione per sito, parametro d'ordine
      magn_temp = 0
      do j=1,lenght
		  if (spin_states(i,j) == 1)  magn_temp = magn_temp + 1
          if (spin_states(i,j) == 0)  magn_temp = magn_temp - 1
      end do
  sum_z = sum_z + abs(dble(magn_temp)) * abs(psi(i))**2              !! pseudo magn (=/0, calcolata con bra e ket)                                                             
  end do														  !! l'altra longitudinale fa 0 con condizioni periodiche
  print *,'Symmetry-broken magnetization along Z:   ', sum_z/lenght
  
  magX = 0.d0
  magZ = 0.d0
     do i = 1,dim_hilbert
         sum_x = 0.d0
         do j = 1,lenght
             if (spin_states(i,j) == 1)   magZ(j) = magZ(j) + abs(Psi(i))**2
             if (spin_states(i,j) == 0)   magZ(j) = magZ(j) - abs(Psi(i))**2
      
             jtemp = lenght +1 -j 
             if (spin_states(i,jtemp) == 1)   temp = i -2**(j-1)
             if (spin_states(i,jtemp) == 0)   temp = i +2**(j-1)
             sum_x = sum_x + dconjg(psi(i)) * psi(temp)
         end do
     if (abs(aimag(sum_x)) > 1.d-10)  stop 'Non real magnetization'    !!check per la routine che usa numeri complessi
     magX(i) = dreal(sum_x)
     end do
  
end subroutine total_magn

subroutine davidson(n, Kmax, Eigenvalue, EigenVector)
  use parameters

  IMPLICIT NONE
  LOGICAL           UseGuess,Wanted
  INTEGER           kmax,jmax,jmin,maxstep,method,m,l,maxnmv,order,testspace,j,lwork,istate,ii,n,kmaxuser
  DOUBLE PRECISION  tol,lock,targetEn,Norm,emin,etemp
  DOUBLE PRECISION, DIMENSION(Kmax)             :: EigenValue
  DOUBLE COMPLEX,   DIMENSION(n)                :: EigenVector
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: alpha,beta,tmp,residu
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: eivec,zwork


  !!  INIZIALIZATION OF PARAMETERS  !!
  Useguess = .false.
  KMaxUser = KMax
  targetEn = -5.d0*lenght
  tol = 1.d-9    ! Tolerance of the eigensolutions: $\Vert \beta H_{SB} x - \alpha x \vert$
  maxnmv = 100    ! Maximum number of matvecs in cgstab or gmres (very problem dependent; typically choose in [5-100])
  wanted = .true. ! If wanted=.true. then computes the converged eigenvectors
  order = -1      ! Selection criterion for Ritz values:  0 (nearest to target);  -1 (smallest real part)
  IF(order == 0)  testspace = 3 ! put 3 if a reasonable value for target is known, else take 2
  IF(order /= 0)  testspace = 2

  IF (3*KmaxUser <= 20) jmax=20          ! Maximum size of the search space:
  IF (3*KmaxUser >  20) jmax=3*KmaxUser
  jmin=2*KmaxUser                        ! Minimum size of the search space
  maxstep = 1000                         ! Maximum number of Jacobi-Davidson iterations
  lock = 1.d-12                          ! Tracking parameter
  method = 2                             ! Method for the linear equation solver  1: gmres(m)  2: cgstab(l)
  m = 30                                 ! Maximum dimension of searchspace for gmres(m):
  l= 2                                   ! Degree of gmres-polynomial in bi-cgstab(l):
  IF (method == 1) lwork =  4 +  m  + 5*jmax + 3*KmaxUser  ! Size of workspace
  IF (method == 2) lwork = 10 + 6*l + 5*jmax + 3*KmaxUser  !KmaxUser is used since Kmax = 1 gives problems ...!
  !!  END OF INIZIALIZATION  !!

  ALLOCATE (alpha(jmax), beta(jmax), eivec(n,Kmax))
  Alpha=0.d0
  Beta=0.d0
  EiVec=0.d0
  ALLOCATE (tmp(n), residu(n), zwork(n,lwork))
  tmp=0.d0
  residu=0.d0
  zwork=0.d0

  CALL JDQZ(ALPHA, BETA, EIVEC, wanted, n, targetEn, tol, Kmax, jmax, jmin, method, m, l, maxnmv, maxstep, &
            lock, order, testspace, zwork, lwork, UseGuess )

  !     Computes the norms of the residuals:
  DO j = 1, Kmax
     CALL AMUL  ( n, eivec(1,j), residu )
     CALL ZSCAL ( n, beta(j), residu, 1 )
     CALL BMUL  ( n, eivec(1,j), tmp )
     CALL ZAXPY ( n, -alpha(j), tmp, 1, residu, 1 )
  ENDDO
  DEALLOCATE (zwork,tmp,residu)
  Eigenvalue(1:Kmax) = dReal(alpha(1:Kmax)/beta(1:Kmax))

  !     Calculates the smallest eigenvalue (ground state)
  emin=eigenvalue(1)
  istate = 1
  DO ii=2,Kmax
     IF (eigenvalue(ii) < emin) THEN
        emin=eigenvalue(ii)
        istate = ii
     ENDIF
  ENDDO
  IF (istate /= 1) THEN
     etemp=eigenvalue(1)
     eigenvalue(1)=eigenvalue(istate)
     eigenvalue(istate)=etemp
  ENDIF
  DEALLOCATE (alpha,beta)

!  print *,'istate',istate
!  Chooses the eigenvector corresponding to the selected eigenvalue
  EigenVector = eivec(:,istate)
  Norm = Sum(dConjg(EigenVector)*(EigenVector))
  EigenVector = EigenVector/(Norm**0.5d0)
  DEALLOCATE (eivec)

end subroutine Davidson
 
subroutine amul(n, psiIn, psiOut)     !!ausiliaria a Davidson, calcola |psiOut> = H*|psiIn>
  use parameters
  implicit none 
  integer  :: n, i, j, jtemp, temp_amul
  double complex, dimension(dim_hilbert), intent(in)  :: psiIn
  double complex, dimension(dim_hilbert), intent(out) :: psiOut

  psiOut = (0.d0,0.d0)

  do i=1,dim_hilbert         !!interazione -sigma_z*sigma_z
	  do j=1,lenght-1
		  if (spin_states(i,j) == spin_states(i,j+1))    psiOut(i) = psiOut(i) - psiIn(i)
		  if (spin_states(i,j) /= spin_states(i,j+1))    psiOut(i) = psiOut(i) + psiIn(i)
      end do
  end do
  
  if (periodic) then  			 !!condizioni periodiche
      do i=1,dim_hilbert
          if (spin_states(i,lenght) == spin_states(i,1))    psiOut(i) = psiOut(i) - psiIn(i)
          if (spin_states(i,lenght) /= spin_states(i,1))    psiOut(i) = psiOut(i) + psiIn(i)
      end do
  end if

  do i=1,dim_hilbert      !!campo longitudinale -sigma(z)
       do j=1,lenght
          if (spin_states(i,j) == 1)   psiOut(i) = psiOut(i) - lambdaAmul*psiIn(i)
          if (spin_states(i,j) == 0)   psiOut(i) = psiOut(i) + lambdaAmul*psiIn(i)
       end do
  end do

  do i = 1,dim_hilbert         !!campo longitudinale -sigma(x)
      do j=1,lenght
          jtemp = lenght +1 -j 
          if (spin_states(i,jtemp) == 1)   temp_amul = i -2**(j-1)    !!conversione da binario, ricontrollare
          if (spin_states(i,jtemp) == 0)   temp_amul = i +2**(j-1)
        psiOut(temp_amul) = psiOut(temp_amul) + g*psiIn(i)
      end do
  end do

end subroutine amul	 

subroutine bmul(neq, q, r ) 		  !!ausiliaria a davidson
  implicit none
  integer :: neq
  double complex :: q(neq),r(neq)
    
  r=q

end subroutine bmul

subroutine Guess(N,X)    			  !!due subroutine ausiliarie non implementate sotto, servono per compilare
  implicit none
  integer :: n
  double complex :: X( * )

  X(1:n) = X(1:n)!PsiGuess(1:n) ! NOT IMPLEMENTED

end subroutine Guess

subroutine precon (neq,psi) 
  !     This subroutine computes $\psi = K \psi$, where $K$ is a matrix which mimics the action of 
  !     $(H - \tau \mathbb{1})^{-1}$. Here H is the Hamiltonian to be diagonalized, $\tau$ is the target 
  !     of the Davidson method, namely the value near which the eigenvalues are sought.
  !     A zeroth order approximation is typically used: $K_{i,j} = \delta_{i,j} (H_{i,i}-\tau)^{-1}$
  implicit none
  integer :: neq
  double complex :: psi(neq)

  psi=psi

end subroutine precon 
		

