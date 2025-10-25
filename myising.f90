program ising_quantum  		!! test per modulo 2 metodi, ground state con davidson
     implicit none 
     
     logical :: periodic
     integer :: lenght, dim_hilbert
     integer :: i, j, itemp
     integer, parameter :: kmax = 3          !!parametro davidson, da guardare
     integer, dimension(:,:), allocatable :: spin_states
     real(8) :: g, lambda, lambdaAmul, cputime
     real(8), dimension(:), allocatable :: en, magn_long, magn_trasv      !!long=z, trasv=x
     double complex, dimension(:), allocatable :: ground_state     
    
     open(unit=10, file='parameters_in.dat', status='old')     !!leggo i parametri rilevanti
		read (10,*) lenght 		!!lunghezza catena
		read (10,*) g 		!!campo trasverso
		read (10,*) lambda 		!!campo longitudinale
		read (10,*) periodic    !!condizioni al bordo periodiche se .true.
	 close(10)
	 
	 dim_hilbert = 2**lenght				!!dimensione dello spazio di hilbert(= dim base computazionale)
	 allocate(spin_states(dim_hilbert, lenght))
	 spin_states = 0
	 
	 do i=1,dim_hilbert             !!inizializzo gli spin e li codifico in binario
		 itemp = i-1                !!da 0 a 7 in binario, codifica al contrario del goat
		 do j=1,lenght
			 spin_states(i,L+1-j)=mod(itemp, 2)
			 itemp = itemp/2
		 end do
	 !print*, i,  spin_states(i, :)      !!verifica stati
	 end do
		 
	 allocate(ground_state(dim_hilbert), magn_long(lenght), magn_trasv(lenght))
	 allocate(en(kmax))
	 lambdaAmul = lambda
	 
	 call(davidson(dim_hilbert, kmax, en, ground_state))    !!magia: trovo il ground state (energia e funzione d'onda)
	 call total_magn(ground_state, magn_trasv, magn_long)   !!calcolo di magnetizzazione per il ground state 
	 
	 print *,'First three energy levels:   ', en(1:3)
	 
	 
	 call cpu_time(cputime)
     print *,'CPU time for the full process:', cputime
	 
	 stop
contains 

subroutine Davidson(n, Kmax, Eigenvalue, EigenVector)
  implicit none
  logical           UseGuess,Wanted
  integer           kmax,jmax,jmin,maxstep,method,m,l,maxnmv,order,testspace,j,lwork,istate,ii,n,kmaxuser
  double precision  tol,lock,targetEn,Norm,emin,etemp
  double precision, dimension(Kmax)             :: EigenValue
  double complex,   dimension(n)                :: EigenVector
  double complex,   dimension(:),   allocatable :: alpha,beta,tmp,residu
  double complex,   dimension(:,:), allocatable :: eivec,zwork


  !!  INIZIALIZATION OF PARAMETERS  !!
  Useguess = .false.
  KMaxUser = KMax
  targetEn = -5.d0*lenght
  tol = 1.d-9    ! Tolerance of the eigensolutions: $\Vert \beta H_{SB} x - \alpha x \vert$
  maxnmv = 100    ! Maximum number of matvecs in cgstab or gmres (very problem dependent; typically choose in [5-100])
  wanted = .true. ! If wanted=.true. then computes the converged eigenvectors
  order = -1      ! Selection criterion for Ritz values:  0 (nearest to target);  -1 (smallest real part)
  if(order == 0)  testspace = 3 ! put 3 if a reasonable value for target is known, else take 2
  if(order /= 0)  testspace = 2

  if (3*KmaxUser <= 20) jmax=20          ! Maximum size of the search space:
  if (3*KmaxUser >  20) jmax=3*KmaxUser
  jmin=2*KmaxUser                        ! Minimum size of the search space
  maxstep = 1000                         ! Maximum number of Jacobi-Davidson iterations
  lock = 1.d-12                          ! Tracking parameter
  method = 2                             ! Method for the linear equation solver  1: gmres(m)  2: cgstab(l)
  m = 30                                 ! Maximum dimension of searchspace for gmres(m):
  l= 2                                   ! Degree of gmres-polynomial in bi-cgstab(l):
  if (method == 1) lwork =  4 +  m  + 5*jmax + 3*KmaxUser  ! Size of workspace
  if (method == 2) lwork = 10 + 6*l + 5*jmax + 3*KmaxUser  !KmaxUser is used since Kmax = 1 gives problems ...!
  !!  END OF INIZIALIZATION  !!

  allocate (alpha(jmax), beta(jmax), eivec(n,Kmax))
  Alpha=0.d0
  Beta=0.d0
  EiVec=0.d0
  allocate (tmp(n), residu(n), zwork(n,lwork))
  tmp=0.d0
  residu=0.d0
  zwork=0.d0

  call JDQZ(ALPHA, BETA, EIVEC, wanted, n, targetEn, tol, Kmax, jmax, jmin, method, m, l, maxnmv, maxstep, &
            lock, order, testspace, zwork, lwork, UseGuess )

  !     Computes the norms of the residuals:
  do j = 1, Kmax
     call AMUL  ( n, eivec(1,j), residu )
     call ZSCAL ( n, beta(j), residu, 1 )
     call BMUL  ( n, eivec(1,j), tmp )
     call ZAXPY ( n, -alpha(j), tmp, 1, residu, 1 )
  end do
  deallocate (zwork,tmp,residu)
  Eigenvalue(1:Kmax) = dReal(alpha(1:Kmax)/beta(1:Kmax))

  !     Calculates the smallest eigenvalue (ground state)
  emin=eigenvalue(1)
  istate = 1
  do ii=2,Kmax
     if (eigenvalue(ii) < emin) then
        emin=eigenvalue(ii)
        istate = ii
     end if
  end do
  if (istate /= 1) then
     etemp=eigenvalue(1)
     eigenvalue(1)=eigenvalue(istate)
     eigenvalue(istate)=etemp
  end if
  deallocate (alpha,beta)

!  print *,'istate',istate
!  Chooses the eigenvector corresponding to the selected eigenvalue
  EigenVector = eivec(:,istate)
  Norm = Sum(dConjg(EigenVector)*(EigenVector))
  EigenVector = EigenVector/(Norm**0.5d0)
  deallocate (eivec)

end subroutine Davidson 
	 
subroutine AMUL(n, psiIn, psiOut)               !!ausiliaria a Davidson, calcola |psiOut> = H*|psiIn>
  implicit none 
  integer  :: n, i, j, temp
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
          if (spin_states(i,lenght) /= spin_states(i,1))    psiOut(ii) = psiOut(ii) + psiIn(i)
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
          if (spin_states(i,j) == 1)   temp = i -2**(j-1)
          if (spin_states(i,j) == 0)   temp = i +2**(j-1)
        psiOut(temp) = psiOut(temp) + g*psiIn(i)
      end do
  end do

end subroutine AMUL	 

subroutine BMUL(neq, q, r ) !  ausiliaria a davidson
  implicit none
  integer :: neq
  double complex :: q(neq),r(neq)
    
  r=q

end subroutine BMUL

subroutine total_magn(psi, magX, magZ)
  implicit none 
  real(8), dimension(lenght) :: magX, magZ
  real(8) :: sum_z
  double complex, dimension(dim_hilbert) :: psi
  double complex :: sum_x
  integer :: i, j, temp, magn_temp
  
  sum_z = 0                  !!magnetizzazione per sito, parametro d'ordine
  do i=1, dim_hilbert
      magn_temp = 0
      do j=1,lenght
		  if (spin_states(i,j) == 1)  magn_temp = magn_temp + 1
          if (spin_states(i,j) == 0)  magn_temp = magn_temp - 1
      end do
  sum_z = sum_z + abs(dble(magn_temp))*abs(ground_state(i))**2     !!? what
  end do
  print *,'Symmetry-broken magnetization along Z:   ', sum_z/lenght
  
end subroutine total_magn
  

end program ising_quantum
		 
	 
	 
	 
		

