module constants
    implicit none
    integer, parameter::n=20,simu=2000,dim=3,ewma=50
    real,parameter:: ZIPpro_ocewma=0.10,ZIPpro_occo=0.90,pi=0.90
end module constants
program main
    use IMSL
    use constants
    implicit none
    integer i,j
    real limit,arl,std
    real ZIPpro_oc(dim,dim),ZIPpro(dim,dim),IC_ZIPpro_oc(dim,dim),miu_exp(dim),Inverse_Cov_exp(dim,dim),shift(dim,dim)
    open(10,file="D:/Multivariate_Poisson/Simulations/MN-GLRT/ZIPpro-1.0-2000-30.txt")

    call parameter_setting(ZIPpro)
    call limit_search(ZIPpro,Limit,arl,std)
    !limit=   1.140625        
	!arl= 197.4545        
    write(10,*),"Limit is",Limit
	write(10,*),"arl is",arl
    shift=0.0
    do i = 1, 9, 1
         if ( i==5 ) then
           cycle 
        end if
        write(10,*),"shifts in (1,1)"
        shift(1,1)=(ZIPpro(1,1)+ZIPpro(1,2)+ZIPpro(1,3))*0.1*(i-5)
        call OCARL(ZIPpro,Limit,shift,arl,std)
        write(10,*),"shift is",0.1*(i-5)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do

    shift=0.0
    do i = 1, 9, 1
         if ( i==5 ) then
           cycle 
        end if
        write(10,*),"shifts in (1,2)"
        shift(1,2)=ZIPpro(1,2)*0.2*(i-5)
        shift(2,1)=shift(1,2)
        shift(1,1)=-shift(1,2)
        shift(2,2)=-shift(1,2)
        call OCARL(ZIPpro,Limit,shift,arl,std)
        write(10,*),"shift is",0.2*(i-5)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do


    shift=0.0
    do i = 1, 9, 1
         if ( i==5 ) then
           cycle 
        end if
        write(10,*),"shifts in (1,3)"
        shift(1,3)=ZIPpro(1,3)*0.2*(i-5)
        shift(3,1)=shift(1,3)
        shift(1,1)=-shift(1,3)
        shift(3,3)=-shift(1,3)
        call OCARL(ZIPpro,Limit,shift,arl,std)
        write(10,*),"shift is",0.2*(i-5)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do


    stop
end program main
subroutine parameter_setting(lamda_0)
    use constants
    use IMSL
    implicit none
    real lamda_0(dim,dim)
    integer i,j

    do i = 1, dim, 1
        do j = 1, dim, 1
            lamda_0(i,j)=1.0/(abs(i-j)+1)
        end do       
    end do 
    return 
    
end subroutine parameter_setting


subroutine sample_generation(lamda_0,sample)
    use IMSL
    use constants
    implicit none
    real lamda_0(dim,dim)
    integer sample(dim,n),samplegeneration(dim,dim,n)
    integer i,j,k

    samplegeneration=0
    do i = 1, dim, 1
            do j = i, dim, 1
                if ( pi<1.0 ) then
                    call RNBIN(n,1,pi,samplegeneration(i,j,:))
                else
                    samplegeneration(i,j,:)=1
                end if
                do k = 1, n, 1
                    if ( samplegeneration(i,j,k)==1 ) then
                        call RNPOI(1,lamda_0(i,j),samplegeneration(i,j,k))
                    end if   
                end do
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
end subroutine sample_generation

subroutine MM_Estimation(sample,miu_ewma,Cov_ewma)
    use constants
    implicit none
    integer i,j,d
    integer sample(dim,n)
    real miu_ewma(dim),Cov_ewma(dim,dim),Cov(dim,dim)

    do i = 1, dim, 1
        miu_ewma(i)=ZIPpro_occo*miu_ewma(i)+ZIPpro_ocewma*sum(sample(i,:))/n
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
    Cov_ewma=ZIPpro_ocewma*Cov+ZIPpro_occo*Cov_ewma 
    return
end subroutine MM_Estimation


subroutine Statistic_meanandcovII(miu_test,miu_exp,Cov_test,Inverse_Cov_exp,Test)
        use IMSL
        use constants
        implicit none
        real Cov_test(dim,dim)
        real Inverse_Cov_exp(dim,dim)
        real Tem_matrix(dim,dim),identy(dim,dim)
        real Test,Trace
        integer i
        real miu_test(1,dim),miu_exp(1,dim),Test_matrix(1,1)


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



subroutine OCARL(ZIPpro,Limit,shift,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real ZIPpro_oc(dim,dim),ZIPpro(dim,dim),shift(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu)
    real Limit,Test,arl,st,std
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 
    !call parameter_settingI(ZIPpro)
    !write(10,*),"IC parameter is",ZIPpro
    
    ZIPpro_oc=ZIPpro+shift
    miu_exp=0.0
    Cov_exp=pi*ZIPpro   
    do i = 1, dim, 1
        do j = 1, dim, 1
        miu_exp(i)=miu_exp(i)+Cov_exp(i,j)         
        end do
    end do

    do i = 1, dim, 1
        do j = 1, dim, 1
        Cov_exp(i,j)=Cov_exp(i,j)*(1.0+(1.0-pi)*ZIPpro(i,j))         
        end do
    end do

    do i = 1, dim, 1
    do j = 1, dim, 1
        if (j/=i)then
        Cov_exp(i,i)=Cov_exp(i,i)+Cov_exp(i,j)
        end if          
        end do
    end do
    Inverse_Cov_exp= .i. Cov_exp

    i=1
    do while(i<=simu)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call sample_generation(ZIPpro,sample)
    call MM_Estimation(sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcovII(miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
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
    !print*,ZIPpro
    !print*,"OK"
    do while (Test<Limit)
    call sample_generation(ZIPpro_oc,sample)
    call MM_Estimation(sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcovII(miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
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

subroutine limit_search(ZIPpro,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j
  integer sample(dim,n)
    real ZIPpro(dim,dim),shift(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu)
    real Limit,Test,arl,st,std,limitl,limitR
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 
    !call parameter_settingI(ZIPpro)
    !write(10,*),"IC parameter is",ZIPpro
    limitl=1.3
    limitR=1.7
    miu_exp=0.0
    Cov_exp=pi*ZIPpro   
    do i = 1, dim, 1
        do j = 1, dim, 1
        miu_exp(i)=miu_exp(i)+Cov_exp(i,j)         
        end do
    end do

    do i = 1, dim, 1
        do j = 1, dim, 1
        Cov_exp(i,j)=Cov_exp(i,j)*(1.0+(1.0-pi)*ZIPpro(i,j))         
        end do
    end do

    do i = 1, dim, 1
    do j = 1, dim, 1
        if (j/=i)then
        Cov_exp(i,i)=Cov_exp(i,i)+Cov_exp(i,j)
        end if          
        end do
    end do
    !print*,Cov_exp
    !print*,miu_exp
    !stop
    Inverse_Cov_exp= .i. Cov_exp
    do while (abs(arl-200.0) >= 4.0)
    Limit = (LimitL + LimitR) / 2.0
    i=1
    do while(i<=simu)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call sample_generation(ZIPpro,sample)
    call MM_Estimation(sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcovII(miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
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
    call sample_generation(ZIPpro,sample)
    call MM_Estimation(sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcovII(miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
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



