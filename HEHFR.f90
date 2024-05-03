program HEHFR
  implicit none
  INTEGER :: i, j, k, l, iter, m
  INTEGER, PARAMETER :: n= 2 !!Número atómico para iterar
  DOUBLE PRECISION, PARAMETER :: Z=2.0D0, re=-2.904D0, c=137.03599D0 !! Número atómico
  DOUBLE PRECISION, dimension(n) :: zeta !! Coeficiente de Slater
  DOUBLE PRECISION, dimension(n,n) :: PL, CLS, CSL, PS, G
  DOUBLE PRECISION, dimension(2*n,2*n) :: S, nc, F, FP, CP, X, X_dag, tem, P, oldP
	DOUBLE PRECISION, dimension(n,n,n,n) :: T  !!Matriz de 4 dimensiones
	DOUBLE PRECISION :: su, mul, E, Delta, crit, er!! Valores
 	
  zeta(1)=1.45363D0 !!Parametro de slater 1
  zeta(2)=2.91093D0 !!Parametro de slater 2

  DO i= 1, n
    DO j= 1, n
			P(i,j)=0.0D0 !!Inicialización de la matriz de solapamiento
    END DO
  END DO

	DO i= 1, n
    DO j= 1, n
			S(i,j)=0.0D0 !!Inicialización de la matriz de solapamiento
    END DO
  END DO

  DO i= 1, n
    DO j= 1, n
      S(i,j)=(8.0D0*(zeta(i)*zeta(j))**(1.5D0))/(zeta(i)+zeta(j))**(3.0D0) !!Matriz de ortogonalidad large
    END DO
  END DO
  
  DO i= (n+1), 2*n
    DO j= (n+1), 2*n
      		S(i,j)=(-4.0D0*(zeta(i-n)*zeta(j-n))**(1.5D0))/(zeta(i-n)+zeta(j-n))**(3.0D0) !!Matriz de ortogonalidad small
    END DO
  END DO

  call find_X_and_X_dagger(S, X, X_dag, n) !!Llama la subrutina para calcular la matriz de transformación
  
  
  DO i= 1, n
    DO j= 1, n
      CLS(i,j)= (4*zeta(i)**(2.5D0)*zeta(j)**(1.5D0))/(zeta(i)+zeta(j))**(3.0D0) !!Matriz de energía cinética large-small
    END DO
  END DO
  
  DO i= 1, n
    DO j= 1, n
      CSL(i,j)= (4*zeta(j)**(2.5D0)*zeta(i)**(1.5D0))/(zeta(i)+zeta(j))**(3.0D0) !!Matriz de energía cinética small-large
    END DO
  END DO
  
  DO i= 1, n
    DO j= 1, n
      PL(i,j)= (-4*Z*(zeta(j)*zeta(i))**(1.5D0))/(zeta(i)+zeta(j))**(2.0D0) !!Matriz de energía potencial large
    END DO
  END DO
  
  DO i= 1, n
    DO j= 1, n
      PS(i,j)= (-4.0D0*(zeta(i)*zeta(j))**(1.5D0))/(zeta(i)+zeta(j))**(3.0D0) - &
      (Z*(zeta(i)*zeta(j))**(1.5D0))/((c*(zeta(i)+zeta(j)))**(2.0D0))!!Matriz de energía potencial small
    END DO
  END DO
  
  DO i= 1, n
    DO j= 1, n
      DO k= 1, n
      	DO l= 1, n
      		su= zeta(i)+zeta(j)+zeta(k)+zeta(l)
      		mul= zeta(i)*zeta(j)*zeta(k)*zeta(l)
      		T(i,j,k,l)= 16.0D0*((mul)**(1.5D0))*((2/((zeta(i)+zeta(k))**(3.0D0)*(zeta(j)+zeta(l))**(2.0D0))) - & 
          (2/((zeta(i)+zeta(k))**(2.0D0)*(su)**(3.0D0)))- & 
          (2/((zeta(i)+zeta(k))**(3.0D0)*(su)**(2.0D0)))) !!Matriz de energía de 2 electrón
      	END DO
      END DO
    END DO
  END DO
  	
  
  Open(12, file="energías", status='unknown')
  write(12,*) '*********************************************************'
  write(12,15)
  15 format (4x, 'Iteración', 2x, 'E_total u.a', 2x, 'E_real u.a', 2x, 'Error %',/)
  
  iter=0
  crit= 1.0D-10 !!Criterio de estabilidad
  1 iter=iter+1
  
  DO i= 1, n
    DO j= 1, n
      DO k= 1, n
      	DO l= 1, n
      		G(i,j) = G(i,j) + P(k,l)*(T(i,k,j,l) - 0.5D0*T(i,k,l,j))!!!!Aporte interacción electrón-electrón
      	END DO
      END DO
    END DO
  END DO
  
  do i = 1, 2
        ! Imprime una fila completa en cada iteración
        write(*, '(4F10.3)') (G(i, j), j = 1, 2)
    end do
  
  DO i= 1, n
    DO j= 1, n
      	F(i,j)=PL(i,j) + G(i,j)!!Matriz de Fock Large
    END DO
  END DO
  
  DO i= (n+1), 2*n
    DO j= (n+1), 2*n
      F(i,j)=PS(i-n,j-n)!!Matriz de Fock Small
    END DO
  END DO
  
  DO i= 1, n
    DO j= (n+1), 2*n
      F(i,j)=CLS(i,j-n)!!Matriz de Fock Large-Small
    END DO
  END DO
  
  DO i= (n+1), 2*n
    DO j= 1, n
      F(i,j)=CSL(i-n,j)!!Matriz de Fock Small-large
    END DO
  END DO
  
	E=0.0D0
  DO i= 1, n
    DO j= 1, n
      	E=E+ P(i,j)*PL(i,j)!!Matriz de Fock Large
    END DO
  END DO
  DO i= 1, n
    DO j= 1, n
      	E=E+ 0.5d0*P(i,j)*G(i,j)!!Matriz de Fock Large
    END DO
  END DO
  DO i= (n+1), 2*n
    DO j= (n+1), 2*n
      E=E+P(i,j)*PS(i-n,j-n)!!Matriz de Fock Small
    END DO
  END DO
  DO i= 1, n
    DO j= (n+1), 2*n
      E=E+P(i,j)*CLS(i,j-n)!!Matriz de Fock Large-Small
    END DO
  END DO
  DO i= 1, n
    DO j= (n+1), 2*n
      E=E-P(i,j)*CSL(i,j-n)!!Matriz de Fock Large-Small
    END DO
  END DO
  
	er=(abs(E+2.904D0)/2.904D0)*100 !!Error relativo porcentual
  
  
	
	tem=matmul(X_dag,F)
	FP=matmul(tem,X) !!Transformación de F a la base ortonormal
	
	call autovector(FP, CP, n) !!Calculo de los autovectores de F en el sistema primado 
	
	nc=0.0D0
	nc=matmul(X,transpose(CP)) !!Nuevos coeficientes de peso
	
	oldP=P
	
	DO i= 1, n
    DO j= 1, n
    	DO k=1, n
		    P(i,j)=(2.0D0*nc(k,i)*nc(j,k)) !!Densidad de probabilidad Large
		  END DO
		END DO
	END DO
	
  
  Delta=0.0D0
  DO i= 1, n
    DO j= 1, n
      Delta = Delta + ABS(oldP(i,j)-P(i,j)) !!Diferencia poblacional
    END DO
  END DO
	
  write (12,18) iter, E, re, er
  18 format(4x,I3,8x,F9.6,3x,F9.6,3x,F9.6)
  
  
  if (Delta .gt. crit) then
		go to 1 !!Retorna el ciclo
	end if
	return
end program HEHFR

subroutine autovector(FP, CP, n) !!Subrutina para calcula autovectores
   	use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    !!Parametros para llamar la función dsyevd
    integer, intent(in) :: n
    real(dp), intent(in) :: FP(2*n, 2*n)
    real(dp), intent(out) :: CP(2*n, 2*n)
    real(dp), dimension(:), allocatable :: work
    real(dp), dimension(2*n) :: eigenvalues
    integer, dimension(:), allocatable :: iwork
    integer :: info, lwork, liwork
    integer :: i,m,j
		
		m=2*n 
    ! Calculamos los tamaños de los arrays de trabajo
    lwork = -1
    liwork = -1
    allocate(work(1))
    allocate(iwork(1))
		
		 ! Llamada de trabajo para dsyevd para calcular el tamaño óptimo de lwork y liwork
    call dsyevd('V', 'U', m, FP, m, eigenvalues, work, lwork, iwork, liwork, info)

    ! Asignamos los tamaños de los arrays de trabajo
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work)
    deallocate(iwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
		
    ! Llamamos a dsyevd para calcular los valores propios y los vectores propios de S
    call dsyevd('V', 'U', m, FP, m, eigenvalues, work, lwork, iwork, liwork, info)

    if (info /= 0) then
        print *, 'Error en la diagonalizacion: Info =', info
        stop
    end if



    ! Formamos la matriz C en la base ortonormal
    CP = FP

		
    
end subroutine autovector



subroutine find_X_and_X_dagger(S, X, X_dag, n)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: S(2*n, 2*n)
    real(dp), intent(out) :: X(2*n, 2*n), X_dag(2*n, 2*n)
    real(dp), dimension(2*n, 2*n) :: D, SP
    real(dp), dimension(:), allocatable :: eigenvalues
    real(dp), dimension(:,:), allocatable :: eigenvectors
    real(dp), dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: iwork
    integer :: info, lwork, liwork
    integer :: i, j, m

    m = 2*n
    ! Asignar tamaño a los arrays de salida
    allocate(eigenvalues(m))
    allocate(eigenvectors(m, m))

    ! Consulta inicial para determinar los tamaños óptimos de los arrays de trabajo
    lwork = -1
    liwork = -1
    allocate(work(1))
    allocate(iwork(1))
    call dsyevd('V', 'U', m, S, m, eigenvalues, work, lwork, iwork, liwork, info)

    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work)
    deallocate(iwork)
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! Llamada final para calcular los autovalores y autovectores
    call dsyevd('V', 'U', m, S, m, eigenvalues, work, lwork, iwork, liwork, info)
    if (info /= 0) then
        print *, 'Error en dsyevd, info =', info
        stop
    endif
		
    ! Inicializamos la matriz D como cero y asignamos la inversa de la raíz cuadrada
    ! de los valores propios en su diagonal
    D = 0.0_dp
    do i = 1, 2*n
        D(i, i) = 1.0_dp / sqrt(abs(eigenvalues(i)))
    end do
  	
  	
    ! Formamos la matriz X multiplicando los vectores propios de S por D
    
    X = matmul(S, D)
		X = transpose (X)
    ! X_dagger es simplemente la transpuesta de X
    X_dag = transpose(X)
    
    
    ! Liberar memoria
    deallocate(work)
    deallocate(iwork)
end subroutine find_X_and_X_dagger

