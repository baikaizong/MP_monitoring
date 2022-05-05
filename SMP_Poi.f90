module constants
    implicit none
    integer, parameter::n=20,simu=10000,dim=3,ewma=50
    real,parameter:: lamdaewma=0.10,lamdaco=0.90
end module constants
program main
    use IMSL
    use constants
    implicit none
    integer i,j
    real limit,arl,std
    real lamda(dim,dim),lamda_0(dim,dim),IC_lamda(dim,dim),miu_exp(dim),Inverse_Cov_exp(dim,dim),shift(dim,dim)
    open(10,file="F:/Fortran/MP_Fortran/Results/MN_LRT/lamda_0.10-p3-2.0-N20.txt")

    call parameter_settingI(lamda_0)
    call limit_search(lamda_0,Limit,arl,std)
    !limit=  1.312500  

    shift=0.0
    do i = 1, 11, 1
         if ( i==6 ) then
           cycle 
        end if
        write(10,*),"shifts in (1,1)"
        shift(1,1)=lamda_0(1,1)*0.05*(i-6)
        call OCARL(lamda_0,Limit,shift,arl,std)
        write(10,*),"shift is",0.05*(i-6)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do

    shift=0.0
    do i = 1, 11, 1
         if ( i==6 ) then
           cycle 
        end if
        write(10,*),"shifts in (1,2)"
        shift(1,2)=lamda_0(1,2)*0.2*(i-6)
        shift(2,1)=shift(1,2)
        shift(1,1)=-shift(1,2)
        shift(2,2)=-shift(1,2)
        call OCARL(lamda_0,Limit,shift,arl,std)
        write(10,*),"shift is",0.2*(i-6)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do


    shift=0.0
    do i = 1, 11, 1
         if ( i==6 ) then
           cycle 
        end if
        write(10,*),"shifts in (1,3)"
        shift(1,3)=lamda_0(1,3)*0.2*(i-6)
        shift(3,1)=shift(1,3)
        shift(1,1)=-shift(1,3)
        shift(3,3)=-shift(1,3)
        call OCARL(lamda_0,Limit,shift,arl,std)
        write(10,*),"shift is",0.2*(i-6)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do


    stop
end program main
subroutine parameter_settingI(lamda_0)
    use constants
    use IMSL
    implicit none
    real lamda_0(dim,dim)
    integer i,j

    do i = 1, dim, 1
        do j = 1, dim, 1
            lamda_0(i,j)=2.0/(abs(i-j)+1)
        end do       
    end do 
    return 
    return 
    
end subroutine parameter_settingI



subroutine Poisample_generation(lamda_0,sample)
    use IMSL
    use constants
    implicit none
    real lamda_0(dim,dim)
    integer sample(dim,n),samplegeneration(dim,dim,n)
    integer i,j

    samplegeneration=0
    do i = 1, dim, 1
            do j = i, dim, 1
                if ( lamda_0(i,j)/=0.0 ) then
                    call RNPOI(n,lamda_0(i,j),samplegeneration(i,j,:))
                end if
                samplegeneration(j,i,:)=samplegeneration(i,j,:)          
            end do        
    end do

    sample=0.0
    do i = 1, dim, 1
        do j = 1, dim, 1
        sample(i,:)= sample(i,:)+samplegeneration(i,j,:)          
        end do
    end do

    return  
end subroutine Poisample_generation


subroutine MM_Estimation(sample,Cov_ewma)
    use constants
    implicit none
    integer i,j,d
    integer sample(dim,n)
    real miu_ewma(dim),Cov_ewma(dim,dim),Cov(dim,dim)

    do i = 1, dim, 1
        miu_ewma(i)=lamdaco*Cov_ewma(i,i)+lamdaewma*sum(sample(i,:))/n
    end do
    Cov=0.0
    do i = 1, dim, 1
    do j=1,dim,1
        do d = 1, n, 1
            Cov(i,j)=Cov(i,j)+1.0*(sample(i,d)-miu_ewma(i))*(sample(j,d)-miu_ewma(j))
        end do
        Cov(i,j)=Cov(i,j)/n
    end do
    end do
    Cov_ewma=lamdaewma*Cov+lamdaco*Cov_ewma 
    do i = 1, dim, 1
        Cov_ewma(i,i)=miu_ewma(i)
    end do
end subroutine MM_Estimation


subroutine Statistic_meanandcovI(miu_exp,Cov_test,Inverse_Cov_exp,Test)
        use IMSL
        use constants
        implicit none
        real Cov_test(dim,dim)
        real Inverse_Cov_exp(dim,dim)
        real Tem_matrix(dim,dim)
        real miu_exp(dim),miu_test(dim),miu_matrix(1,dim)
        real Test,Trace,tem_var,Test_matrix(1,1)
        integer i
        do i = 1, dim, 1
          miu_test(i)=Cov_test(i,i)
        end do
        miu_matrix(1,:)=miu_test-miu_exp
        Test=0.0
        Tem_matrix=Cov_test .x. Inverse_Cov_exp
        tem_var=det(Tem_matrix)
        Test_matrix=miu_matrix .x. Inverse_Cov_exp .xt. miu_matrix
        Test=1.0*n*(Trace(Tem_matrix)-log(tem_var)-dim+Test_matrix(1,1))
		!print*,Trace(Tem_matrix)-log(tem_var)-dim
		!print*,Test_matrix(1,1)
	!	print*,"---------------------------"
        return
end subroutine Statistic_meanandcovI
subroutine Statistic_meanandcovII(miu_exp,Cov_test,Inverse_Cov_exp,Test)
        use IMSL
        use constants
        implicit none
        real Cov_test(dim,dim)
        real Inverse_Cov_exp(dim,dim)
        real Tem_matrix(dim,dim),identy(dim,dim)
        real Test,Trace
        integer i
        real miu_test(1,dim),miu_exp(1,dim),Test_matrix(1,1)

        do i = 1, dim, 1
          miu_test(1,i)=Cov_test(i,i)
        end do
        Test_matrix=0.0
        Test_matrix=(miu_test-miu_exp) .x. (Inverse_Cov_exp) .xt.  (miu_test-miu_exp)

        identy=0.0
        do i = 1, dim, 1
            identy(i,i)=1.0
        end do

        
        Tem_matrix=Cov_test .x. Inverse_Cov_exp
        Tem_matrix=Tem_matrix-identy
        Tem_matrix=Tem_matrix .x. Tem_matrix
        Test=Trace(Tem_matrix)
        Test=1.0*n*(0.5*Test+Test_matrix(1,1))
        return
end subroutine Statistic_meanandcovII

function Trace(matrixcu)
        use constants
        implicit none
        real matrixcu(dim,dim),trace
        integer i
        Trace=0.0
        do i = 1, dim, 1
            Trace=Trace+matrixcu(i,i)
        end do
        return
end function Trace



subroutine OCARL(lamda_0,Limit,shift,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real lamda(dim,dim),lamda_0(dim,dim),shift(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu)
    real Limit,Test,arl,st,std
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 
    !call parameter_settingI(lamda_0)
    !write(10,*),"IC parameter is",lamda_0
    
    lamda=lamda_0+shift
    Cov_exp=lamda_0
    do i = 1, dim, 1
        do j = 1, dim, 1
        if (j/=i)then
        Cov_exp(i,i)=Cov_exp(i,i)+Cov_exp(i,j)
        end if          
        end do
        miu_exp(i)=Cov_exp(i,i)
    end do
    Inverse_Cov_exp= .i. Cov_exp
    i=1
    do while(i<simu)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call Poisample_generation(lamda_0,sample)
    call MM_Estimation(sample,Cov_ewma)
    call Statistic_meanandcovII(miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
    !print*,Test
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        cycle
    endif   
    !print*,Test
  !stop
    simuarl(i)=0.0
    !print*,"OK"
    do while (Test<Limit)
    call Poisample_generation(lamda,sample)
    call MM_Estimation(sample,Cov_ewma)
    call Statistic_meanandcovII(miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
      !print*,Test
    simuarl(i)=simuarl(i)+1.0
    end do
    print*,i
    print*,simuarl(i)
  !stop
    print*,"****************"
    i=i+1
    end do
    arl=(sum(simuarl))/(1.0*(simu))

   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu-1.0))
  std=std/(sqrt(dble(simu)-1)) 
  return 
end subroutine OCARL

subroutine limit_search(lamda_0,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j
  integer sample(dim,n)
    real lamda_0(dim,dim),shift(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu)
    real Limit,Test,arl,st,std,limitl,limitR
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 
    !call parameter_settingI(lamda_0)
    !write(10,*),"IC parameter is",lamda_0
    limitl=1.0
    limitR=1.2
    Cov_exp=lamda_0
    do i = 1, dim, 1
        do j = 1, dim, 1
        if (j/=i)then
        Cov_exp(i,i)=Cov_exp(i,i)+Cov_exp(i,j)
        end if          
        end do
        miu_exp(i)=Cov_exp(i,i)
    end do
    Inverse_Cov_exp= .i. Cov_exp
    do while (abs(arl-200.0) >= 4.0)
    Limit = (LimitL + LimitR) / 2.0
    i=1
    do while(i<simu)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call Poisample_generation(lamda_0,sample)
    call MM_Estimation(sample,Cov_ewma)
    call Statistic_meanandcovII(miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
    !print*,Test
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        cycle
    endif   
    !print*,Test
  !stop
    simuarl(i)=0.0
    !print*,"OK"
    do while (Test<Limit)
    call Poisample_generation(lamda_0,sample)
    call MM_Estimation(sample,Cov_ewma)
    call Statistic_meanandcovII(miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
      !print*,Test
    simuarl(i)=simuarl(i)+1.0
    end do
    print*,i
    print*,simuarl(i)
  !stop
    print*,"****************"
    i=i+1
    end do
    arl=(sum(simuarl))/(1.0*(simu))
    write(10,*),"limit is",limit
    write(10,*),"arl is",arl
    if (arl < 200.0) then
      LimitL = Limit
     else
      LimitR = Limit
     end if
    enddo
   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu-1.0))
  std=std/(sqrt(dble(simu)-1)) 
  return 
end subroutine limit_search



