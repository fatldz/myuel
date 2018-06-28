c
c simple 2-d linear beam element with generalized section properties
c
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predf,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)
c
       include 'aba_param.inc'
c
       dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),svars(*),props(*),
     1 energy(7),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflag(4),jprops(*)
c
       dimension b(2,7),gauss(2)
c     b矩阵和高斯点矩阵
       parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     1 six=6.d0,eight=8.d0,twelve=12.d0)
       data gauss/.211324865d0,.788675135d0/
c
c calculate length and directional cosines
c
       dx=coords(1,2)-coords(1,1)
       dy=coords(2,2)-coords(2,1)
       dl2=dx**2+dy**2
       dl=sqrt(dl2)
c dl是单元长度
       hdl=dl/two
       acos=dx/dl
       asin=dy/dl
c 计算方向余弦
c
c
       nsvint=nsvars/2
c 每个积分点上的状态变量个数
       do k1=1,7
	   rhs(k1,1)=zero
	   do k2=1,7
	     amatrx(k1,k2)=zero
	   end do
       end do
c initialize rhs and lhs
c
c loop over intergration points
c 以下对高斯积分点循环
       do kintk=1,2
         g=gauss(kintk)
c make b-matrix
c
         b(1,1)=(-three+four*g)*acos/dl
	   b(1,2)=(-three+four*g)*asin/dl
	   b(1,3)=zero
	   b(1,4)=(-one+four*g)*acos/dl
	   b(1,5)=(-one+four*g)*asin/dl
	   b(1,6)=zero
	   b(1,7)=(four-eight*g)/dl
	   b(2,1)=(-six+twelve*g)*(-asin)/dl2
         b(2,2)=(-six+twelve*g)*acos/dl2
	   b(2,3)=(-four+six*g)/dl
	   b(2,4)=(six-twelve*g)*(-asin)/dl2
	   b(2,5)=(six-twelve*g)*acos/dl2
	   b(2,6)=(-two+six*g)/dl
	   b(2,7)=zero     
c
c calculate incremental strains and curvatures
c
         eps=zero
c 应变
         deps=zero
c 增量应变
         cap=zero
c 曲率
      dcap=zero
c 增量曲率      
c 以下语句根据b矩阵和位移（转角）计算应变，增量应变，曲率和增量曲率
      do k=1,7
           eps=eps+b(1,k)*u(k)
           deps=deps+b(1,k)*du(k,1)
           cap=cap+b(2,k)*u(k)
           dcap=dcap+b(2,k)*du(k,1)
      end do 
c
c call constitutive routine ugenb 
c 调用本够关系子程序得到力和弯矩
c
         isvint=1+(kintk-1)*nsvint
c kintk积分点，nsvint状态变量个数
         bn=zero
c 轴力
         bm=zero
c 弯矩
         daxial=zero
c 梁轴向刚度
         dcoupl=zero
c 轴向刚度和弯曲刚度的耦合项
         call ugenb(bn,bm,daxial,dbend,dcoupl,eps,deps,cap,dcap,
     1              svars(isvint),nsvint,props,nprops)
c
c assemble rhs and lhs
c 对RHS和AMARTX赋值
c
         do k1=1,7
           rhs(k1,1)=rhs(k1,1)-hdl*(bn*b(1,k1)+bm*b(2,k1))
c  Bt*D*B*u,rhs就是节点力向量
           bd1=hdl*(daxial*b(1,k1)+dcoupl*b(2,k1))
           bd2=hdl*(dcoupl*b(1,k1)+dbend*b(2,k1))
           do k2=1,7
             amatrx(k1,k2)=amatrx(k1,k2)+bd1*b(1,k2)+bd2*b(2,k2)
c amatrx即是最后的刚度矩阵
           end do
         end do
       end do
c
       return 
       end


      subroutine ugenb(bn,bm,daxial,dbend,dcoupl,eps,deps,cap,dcap,
     1                 svint,nsvint,props,nprops)
c
       include 'aba_param.inc'
c
       parameter(zero=0.d0,twelve=12.d0)
c
       dimension svint(*),props(*)
c
c variables to be defined by the user
c
c bn-axial force
c bm-bending moment
c daxial-current tangent axial stiffness
c dbend-current tangent bending stiffness
c dcoupl-tangent coupling term
c
c variables to be updated
c
c svint-state variables for this intergration point
c
c variables passed in for information
c
c eps-axial strain
c deps-incremental axial strain
c cap-curvature change
c dcap-incremental curvature change 
c props-element properties
c nprops= # element properties 
c nsvint- #state variables
c
c current assumption
c
c props(1)-section height
c props(2)-section width
c props(3)-Young's modulus
c
       h=props(1)
       w=props(2)
       E=props(3)
c
c formulate linear stiffness
c
       daxial=E*h*w
       dbend=E*w*h**3/twelve
       dcoupl=zero
c
c calculate axial force and moment
c  
       bn=svint(1)+daxial*deps
       bm=svint(2)+dbend*dcap
c
c store internal variables 
c
       svint(1)=bn
       svint(2)=bm
       svint(3)=eps
       svint(4)=cap
c
c 每个积分点都有svint(1~4),分别代表轴力，弯矩，轴向应变和曲率
       return
       end
