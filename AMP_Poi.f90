module constants
    implicit none
    integer, parameter::n=20,simu=2000,dim=3,ewma=20,window_size=200
    real,parameter:: lamdaewma=0.10,lamdaco=0.90,epsilon=0.000001
end module constants
program main
    use IMSL
    use constants
    implicit none
    integer i,j
    real*8 limit,arl,std
    real lamda(dim,dim),lamda_0(dim,dim),IC_lamda(dim,dim),miu_exp(dim),Inverse_Cov_exp(dim,dim),shift(dim,dim)
    open(10,file="./reults/AMP_Poi.txt")
    call parameter_setting(lamda_0)

    !call limit_search(lamda_0,Limit,arl,std)
    limit=  0.865000002086163        
	arl= 199.860500000000     
    write(10,*),"Limit is",Limit
	write(10,*),"arl is",arl
    !stop
    shift=0.0
    do i = 1, 11, 1
        if (i==6) then
        cycle
        endif
        write(10,*),"shifts in (1,1)"
        shift(1,1)=lamda_0(1,1)*0.05*(i-6)
        call OCARL(lamda_0,Limit,shift,arl,std)
        write(10,*),"shift is",0.05*(i-6)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do
    shift=0.0
    do i = 1, 17, 1
		 if (i==9) then
		cycle
		endif
        write(10,*),"shifts in (1,2)"
        shift(1,2)=lamda_0(1,2)*0.05*(i-9)
        shift(2,1)=shift(1,2)
        call OCARL(lamda_0,Limit,shift,arl,std)
        write(10,*),"shift is",0.05*(i-9)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
    end do

    shift=0.0
    do i = 1, 17, 1
		 if (i==9) then
		cycle
		endif
        write(10,*),"shifts in (1,3)"
        shift(1,3)=lamda_0(1,3)*0.05*(i-9)
        shift(3,1)=shift(1,3)
        call OCARL(lamda_0,Limit,shift,arl,std)
        write(10,*),"shift is",0.05*(i-9)
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
            lamda_0(i,j)=0.5/(abs(i-j)+1)
        end do       
    end do 
    return 
    
end subroutine parameter_setting



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



subroutine LRT_test(sample,lamda_IC,globaltest)
    use constants
    use IMSL
    implicit none
    integer sample(window_size,dim,n),max_val(dim)
    real lamda_IC(dim,dim),lamda_est(dim,dim),lamda_last(dim,dim)
    real*8 globaltest
    integer i,j,k
    integer iterations
    real*8 FT,likelihoodlast,likelihood,likelihoodnull
    real*8 likelihood_array,S(3),S_ewma(3)
    real*8, allocatable :: pro_table(:,:,:)
    integer tmp(dim,window_size)

    do i = 1, dim, 1
        do j = 1, window_size, 1
            tmp(i,j)=maxval(sample(j,i,:))
        end do
        max_val(i)=maxval(tmp(i,:))
    end do

    allocate(pro_table(max_val(1)+2,max_val(2)+2,max_val(3)+2))
    lamda_est=lamda_IC
    pro_table=0.0
    
    pro_table(2,2,2)=exp(-(lamda_est(1,1)+lamda_est(2,2)+lamda_est(3,3)-lamda_est(1,2)-lamda_est(1,3)-lamda_est(2,3))) 

    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2)= pro_table(2,2,2)*(1.0*((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2)= pro_table(2,2,2)*(1.0*((lamda_est(2,2)-lamda_est(1,2)-lamda_est(2,3))**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i)= pro_table(2,2,2)*(1.0*((lamda_est(3,3)-lamda_est(1,3)-lamda_est(2,3))**(i-2))/FAC(i-2))
    end do 
    !print*,"--------------------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2)=1.0*((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))*pro_table(i-1,j,2)+lamda_est(1,2)*pro_table(i-1,j-1,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j)=1.0*((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))*pro_table(i-1,2,j)+lamda_est(1,3)*pro_table(i-1,2,j-1))/(i-2)
        end do  
    end do 
   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j)=1.0*((lamda_est(2,2)-lamda_est(1,2)-lamda_est(2,3))*pro_table(2,i-1,j)+lamda_est(2,3)*pro_table(2,i-1,j-1))/(i-2)
        end do  
    end do 
    !print*,"--------------------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k)=((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))*pro_table(i-1,j,k)+lamda_est(1,2)*pro_table(i-1,j-1,k)+lamda_est(1,3)*pro_table(i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do

    likelihood=0.0
    do j = 1, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       !print*,pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2)
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
    end do
    if ( j==1 ) then
       likelihood=likelihood_array
    else
       likelihood=lamdaewma*likelihood_array+lamdaco*likelihood 
    end if         
    end do
    likelihoodnull=likelihood 
    S=0.0
    do j = 1, window_size, 1
        do i = 1, 3, 1
           S(i)=1.0*sum(sample(j,i,:))/n         
           if ( j==1 ) then
            lamda_est(i,i)=S(i)
           else
           lamda_est(i,i)=lamdaewma*S(i)+lamdaco*lamda_est(i,i) 
        end if   
        end do
    end do
 
    FT=1.0
    likelihoodlast=-100000.0
    iterations=0
    do while ( FT>epsilon)
    pro_table=0.0
    
    pro_table(2,2,2)=exp(-(lamda_est(1,1)+lamda_est(2,2)+lamda_est(3,3)-lamda_est(1,2)-lamda_est(1,3)-lamda_est(2,3))) 
    !int*,"*********************"   
    !print*,"--------------------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2)= pro_table(2,2,2)*(1.0*((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2)= pro_table(2,2,2)*(1.0*((lamda_est(2,2)-lamda_est(1,2)-lamda_est(2,3))**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i)= pro_table(2,2,2)*(1.0*((lamda_est(3,3)-lamda_est(1,3)-lamda_est(2,3))**(i-2))/FAC(i-2))
    end do 
    !print*,"--------------------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2)=1.0*((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))*pro_table(i-1,j,2)+lamda_est(1,2)*pro_table(i-1,j-1,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j)=1.0*((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))*pro_table(i-1,2,j)+lamda_est(1,3)*pro_table(i-1,2,j-1))/(i-2)
        end do  
    end do 
   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j)=1.0*((lamda_est(2,2)-lamda_est(1,2)-lamda_est(2,3))*pro_table(2,i-1,j)+lamda_est(2,3)*pro_table(2,i-1,j-1))/(i-2)
        end do  
    end do 
    !print*,"--------------------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k)=((lamda_est(1,1)-lamda_est(1,2)-lamda_est(1,3))*pro_table(i-1,j,k)+lamda_est(1,2)*pro_table(i-1,j-1,k)+lamda_est(1,3)*pro_table(i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do

    likelihood=0.0
    !print*,"--------------------"
    !print*,sum(pro_table)
    !print*,"iterations is",iterations
    do j = 1, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       !print*,pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2)
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
    end do
    if ( j==1 ) then
       likelihood=likelihood_array
    else
       likelihood=lamdaewma*likelihood_array+lamdaco*likelihood 
    end if         
    end do
    S_ewma=0.0
    do j = 1, window_size, 1
    S=0.0
    do i = 1, n, 1
         S(1)=S(1)+(lamda_est(1,2)*pro_table(sample(j,1,i)+1,sample(j,2,i)+1,sample(j,3,i)+2)/pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
         S(2)=S(2)+(lamda_est(1,3)*pro_table(sample(j,1,i)+1,sample(j,2,i)+2,sample(j,3,i)+1)/pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
         S(3)=S(3)+(lamda_est(2,3)*pro_table(sample(j,1,i)+2,sample(j,2,i)+1,sample(j,3,i)+1)/pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
    end do
    if ( j==1 ) then
       S_ewma(1)=S(1)
       S_ewma(2)=S(2)
       S_ewma(3)=S(3)
    else
       S_ewma(1)=lamdaewma*S(1)+lamdaco*S_ewma(1)
       S_ewma(2)=lamdaewma*S(2)+lamdaco*S_ewma(2)
       S_ewma(3)=lamdaewma*S(3)+lamdaco*S_ewma(3)
    end if
    enddo
    lamda_last=lamda_est
    !print*,"lamda_est is",lamda_est
    lamda_est(1,2)=S_ewma(1)/n
    lamda_est(2,1)=lamda_est(1,2)
    lamda_est(1,3)=S_ewma(2)/n
    lamda_est(3,1)=lamda_est(1,3)
    lamda_est(2,3)=S_ewma(3)/n
    lamda_est(3,2)=lamda_est(2,3)
    !print*,"likelihood is",likelihood
    FT=abs(likelihood-likelihoodlast)
    if ( likelihood-likelihoodlast<0.0 ) then
      lamda_est = lamda_last
      likelihood=likelihoodlast
      exit
    end if      
    likelihoodlast=likelihood
    iterations=iterations+1
    if ( iterations>100 ) then
        exit
    end if
    end do
    globaltest=2.0*(likelihood-likelihoodnull)
    return
end subroutine LRT_test


subroutine OCARL(lamda_0,Limit,shift,arl,std)
use IMSL
use constants
implicit none
    integer i,j
  integer sample(window_size,dim,n)
    real lamda(dim,dim),lamda_0(dim,dim),shift(dim,dim),lamda_IC(dim,dim)
    real simuarl(simu)
    real*8 Limit,arl,st,std,Test
    !call parameter_setting(lamda_0)
    lamda=lamda_0+shift

    lamda_IC=lamda_0
    do i = 1, dim, 1
        do j = 1, dim, 1
        if (j/=i)then
        lamda_IC(i,i)=lamda_IC(i,i)+lamda_IC(i,j)
        end if          
        end do
    end do

    i=1
    do while(i<simu)
    Test=0.0
    do j = 1, window_size, 1 
    call Poisample_generation(lamda_0,sample(j,:,:))
    end do
    call LRT_test(sample,lamda_IC,Test)
	!print*,Test
    if (Test>Limit) then
        cycle
    endif   
    simuarl(i)=0.0
    do while (Test<Limit)
    sample(1:window_size-1,:,:)=sample(2:window_size,:,:)
    call Poisample_generation(lamda,sample(window_size,:,:))
    call LRT_test(sample,lamda_IC,Test)
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
  integer sample(window_size,dim,n)
    real lamda_0(dim,dim),shift(dim,dim),lamda_IC(dim,dim)
    real simuarl(simu)
    real*8 Limit,arl,st,std,Test,LimitL,LimitR
    !call parameter_setting(lamda_0)
    !write(10,*),"IC parameter is",lamda_0
    lamda_IC=lamda_0
    do i = 1, dim, 1
        do j = 1, dim, 1
        if (j/=i)then
        lamda_IC(i,i)=lamda_IC(i,i)+lamda_IC(i,j)
        end if          
        end do
    end do
    LimitL=0.365
    LimitR=0.370
    do while (abs(arl-200.0) >=4.0)
    Limit = (LimitL + LimitR) / 2.0
    i=1
    do while(i<simu)
    Test=0.0
    do j = 1, window_size, 1 
    call Poisample_generation(lamda_0,sample(j,:,:))
    end do
    call LRT_test(sample,lamda_IC,Test)
	!print*,Test
    if (Test>Limit) then
        cycle
    endif   
    simuarl(i)=0.0
    !print*,"OK"
    do while (Test<Limit)
    sample(1:window_size-1,:,:)=sample(2:window_size,:,:)
    call Poisample_generation(lamda_0,sample(window_size,:,:))
    call LRT_test(sample,lamda_IC,Test)
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