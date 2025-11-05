module parameters			 !!non posso usare contains(subroutines), problemi nel compilare
  logical  :: periodic
  integer  :: lenght, dim_hilbert
  integer, dimension(:,:), allocatable :: spin_states
  real(8) :: g, lambda, lambdaAmul
end module parameters

program energy_gaps			 !!programma per l'analisi dei gap energetici in ising quantistico 1D
	use parameters
	implicit none
	
	integer :: i, j, itemp
	integer :: kmax          !!parametro davidson, autovalore massimo che la routine trova
	real(8) :: gap0, gap1, gap2, cputime
	logical :: save_time, lenght_is_fixed
	real(8), dimension(:), allocatable :: en
	double complex, dimension(:), allocatable :: ground_state
	  
	lenght_is_fixed=.false.              !!se true devo mandare in input il campo, altrimenti la lunghezza    
	save_time = .false.  	 			 !!printa il tempo di calcolo o no 
	periodic = .false.		 		     !!condizioni al bordo
	lambda = 0.d0			 		     !!campo longitudinale (0 per questo progetto)
	
	if (lenght_is_fixed) then            !!inizializzo i parameteri a seconda dello scopo
	    lenght = 10              		 !!lunghezza catena (costante)
	    read *, g 		     	 		 !!campo trasverso (da mandare con script shell)
	else
	    read *, lenght                   !!lunghezza catena (da mandare con script shell)
	 	g = 0.75d0	     	             !!campo trasverso (costante)
    end if  
	
	dim_hilbert = 2**lenght	 !!inizializzo i registri di spin
	allocate(spin_states(dim_hilbert, lenght))
	spin_states = 0
	
	do i=1,dim_hilbert       !!codifica dei registri in binario: in ordine crescente     
		 itemp = i-1         !!con codifica convenzionale da dx
		 do j=1,lenght
			 spin_states(i,lenght+1-j)=mod(itemp, 2)
			 itemp = itemp/2
		 end do
	!print*, i,  spin_states(i, :)      !!verifica stati
	end do
	
	kmax = 4				 
	allocate(ground_state(dim_hilbert))
	allocate(en(kmax))
	
	lambdaAmul = lambda
	call davidson(dim_hilbert, kmax, en, ground_state)    !!trovo il ground state (funzione d'onda)
	if (lenght_is_fixed) then                             !!e i primi kmax livelli energetici
        print *, g, en(1:kmax) 	      
	else 
	    print *, lenght, en(1:kmax)
	end if
	
	if (save_time) then
	    call cpu_time(cputime)
        print *,'Total CPU time:	', cputime
    end if
    
    deallocate(spin_states, ground_state, en)
    
end program energy_gaps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine davidson(n, Kmax, Eigenvalue, EigenVector) !!non toccare
  use parameters

  implicit none
  logical           UseGuess,Wanted
  integer           kmax,jmax,jmin,maxstep,method,m,l,maxnmv,order,testspace,j,lwork,istate,ii,n,kmaxuser
  real(8)           tol,lock,targetEn,Norm,emin,etemp
  real(8), dimension(Kmax)		                :: EigenValue
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

end subroutine davidson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amul(n, psiIn, psiOut)     !!ausiliaria a Davidson, calcola |psiOut> = H*|psiIn>
  use parameters                      !!unica cosa da modificare
  implicit none 
  integer  :: n, i, j, jtemp, temp_amul
  double complex, dimension(dim_hilbert), intent(in)  :: psiIn
  double complex, dimension(dim_hilbert), intent(out) :: psiOut

  psiOut = (0.d0,0.d0)

  do i=1,dim_hilbert           !!interazione -sigma_z*sigma_z
	  do j=1,lenght-1
		  if (spin_states(i,j) == spin_states(i,j+1))    psiOut(i) = psiOut(i) - psiIn(i)
		  if (spin_states(i,j) /= spin_states(i,j+1))    psiOut(i) = psiOut(i) + psiIn(i)
      end do
  end do
  
  if (periodic) then  	       !!condizioni periodiche
      do i=1,dim_hilbert
          if (spin_states(i,lenght) == spin_states(i,1))    psiOut(i) = psiOut(i) - psiIn(i)
          if (spin_states(i,lenght) /= spin_states(i,1))    psiOut(i) = psiOut(i) + psiIn(i)
      end do
  end if

  do i=1,dim_hilbert           !!campo longitudinale -sigma(z)
       do j=1,lenght
          if (spin_states(i,j) == 1)   psiOut(i) = psiOut(i) - lambdaAmul*psiIn(i)
          if (spin_states(i,j) == 0)   psiOut(i) = psiOut(i) + lambdaAmul*psiIn(i)
       end do
  end do

  do i = 1,dim_hilbert         !!campo longitudinale -sigma(x)
      do j=1,lenght
          jtemp = lenght +1 -j 
          if (spin_states(i,jtemp) == 1)   temp_amul = i -2**(j-1)   
          if (spin_states(i,jtemp) == 0)   temp_amul = i +2**(j-1)
        psiOut(temp_amul) = psiOut(temp_amul) + g*psiIn(i)
      end do
  end do

end subroutine amul	 

subroutine bmul(neq, q, r ) 		  !!tre subroutine ausiliarie qua sotto, servono per compilare
  implicit none
  integer :: neq
  double complex :: q(neq),r(neq)
    
  r=q

end subroutine bmul

subroutine Guess(N,X)    			  
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
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
