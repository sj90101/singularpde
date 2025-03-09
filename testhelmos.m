  function testhelmos()
    close all
    format long e
    format compact
    set(0,'DefaultAxesFontSize',12)

    nsub_test = 0; % plot n_sub versus error
    plottiming = 0; % plot timing results and zk versus N_iter
    plotfield = 1; % plot the field and the error
    
%   examples in Bruno & Lintner 2012 Radio Sciences:
%            icase=0 -> line segment
%            icase=1 -> spiral
%   other examples:
%            icase=2 -> corners 
%            icase=3 -> branch points 
%            icase=4 -> maze
    icase = 4
    
%   bctype='d' -> Dirichlet condition
%   bctype='n' -> Neumann condition    
    bctype = 'd';

    helmosopts = [];

    if icase == 2
      helmosopts.ncorner = 1;
      if plotfield == 1
	helmosopts.ncorner = 8;
      end
    end

    if icase == 3
      helmosopts.max_depth = 0;
      if plotfield == 1
	helmosopts.max_depth = 2;
      end
    end
    
    [arclen,zsrc,chnkr,osparams] = compute_curve_length(icase,helmosopts);
    Ltot = sum(arclen);
    ncurve = length(arclen);

    
    if plottiming == 1
      if icase == 1
	tol = 1e-6;
      elseif icase == 2 || icase == 3
	tol = 1e-12;
      end
      Ng = 300;
    elseif plotfield == 1
      tol = 1d-13;
%      Ng = 1000;
      Ng = 3000;
    elseif nsub_test == 1
      tol = eps;
      Ng = 1;
    end

    if plottiming == 1
      if icase == 1
	zk = 100*pi/Ltot % starts from L/lambda = 50
      elseif icase == 2 || icase == 3
	zk = 20*pi/Ltot % starts from L/lambda = 10
      end
    elseif plotfield == 1
      if icase <=3
	zk = 400*pi/Ltot % as in BrunLint12: L/lambda = 200
      end
      if icase == 4
%	zk = 80*pi
%	zk = 40*pi
	zk = 10*pi
      end
    elseif nsub_test == 1
      zk = 3;
    end

    rellen = Ltot*zk/(2*pi);
    disp(['Total length relative to the wavelength = ',num2str(rellen)])


    npanvec = zeros(ncurve,1);
    for i=1:ncurve
      L=arclen(i);
%     number of panels on each curve
%     should be proportional to the length/wavelength
      if tol<1e-14
      	npan0 = round(3*L*abs(zk)/(2*pi));
      elseif tol<1e-11
      	npan0 = round(3*L*abs(zk)/(2*pi));
      elseif tol<1e-5
	if icase == 0
      	  npan0 = round(0.5*L*abs(zk)/(2*pi));
	elseif icase == 1
      	  npan0 = round(L*abs(zk)/(2*pi));
	else
      	  npan0 = round(2.5*L*abs(zk)/(2*pi));
	end
      end
%     number of panels on each curve
      npan0 = max(4,npan0);
      npanvec(i) = npan0;
% should really do some rebalance here to make sure that panels connected
% to the same singular point have roughly the same length?
    end
    npanvec'
    
    if icase == 0
      xylim=[-1.3 1.3 -1.5 1.1];
    elseif icase == 1
      xylim=[-2.5 1.75 -2.75 1.5];
    elseif icase >= 2
      xmin = min(real(zsrc));
      xmax = max(real(zsrc));
      ymin = min(imag(zsrc));
      ymax = max(imag(zsrc));
      lx = xmax-xmin;
      ly = ymax-ymin;

      xcen = (xmax+xmin)/2;
      ycen = (ymax+ymin)/2;

      if icase <= 3
	if lx <= ly
	  ymin = ymin-ly/4;
	  ymax = ymax+ly/4;
	  xmin = xcen-ly*3/4;
	  xmax = xcen+ly*3/4;
	else
	  xmin = xmin-lx/4;
	  xmax = xmax+lx/4;
	  ymin = ycen-lx*3/4;
	  ymax = ycen+lx*3/4;
	end
	xmin = floor(xmin);
	xmax = ceil(xmax);
	ymin = floor(ymin);
	ymax = ceil(ymax);
      end

      if icase == 4
	xmin = floor(xmin)-1;
	xmax = ceil(xmax)+1;
	ymin = floor(ymin)-1;
	ymax = ceil(ymax)+1;
      end
      
%      xylim = [xmin-0.5 xmax ymin ymax+0.5] % icase=2 one corner
%      xylim = [xmin xmax ymin+1 ymax+1] % icase=3 one branch
      xylim = [xmin xmax ymin ymax]
%      xylim=[xmin-lx/4 xmax+lx/4 ymin-ly/4 ymax+ly/4];
    end

%   plot boundary
    plotboundary = 1;
    if plotboundary == 1
      if icase >=3
	cpars=osparams.cpars;
	plot_linesegments(cpars,xylim,icase);
      else
	plot_chunker(chnkr,xylim,icase);
      end
      return
    end

    
%   do nsub test only for a straight line segment    
    if nsub_test == 1
      nsubup = 100;
      if icase == 0
	npanvec = 6;
      end
      
      ztarg=0.17+0.62i;
      if icase == 0
	uref=0.02788626934981090-0.75932847390327920i;
	helmosopts.ztarg = ztarg;
      end
      
      relerr = [];
      for nsub = 1:nsubup
      	nsub
      	osoutput = helmos(icase,tol,bctype,nsub,zk,npanvec,xylim,Ng,helmosopts);
      	unum = osoutput.unum;

	if icase > 0
	  nsub_ref = nsub + 5;
	  npanvec_ref = round(npanvec*1.5);
      	  osoutput_ref = helmos(icase,tol,bctype,nsub_ref,zk,npanvec_ref, ...
				xylim,Ng,helmosopts);
	  uref = osoutput_ref.unum;
	end
	
      	rerr2 = norm(unum-uref)/norm(uref)
      	relerr=[relerr;rerr2];
      end
      nsub_error_plot(relerr,nsubup,nsubup,bctype,icase)
    end

    
%   level of dyadic refinement in the forward recursion for computing R  
     nsub = 40;

    
    if plottiming == 1
      totaltime=[];np=[];gmrestime=[];buildtime=[];evaltime=[];
      relerr=[];rellinferr=[];iter=[];zktot=[];

      nsub_ref = round(nsub*1.5);
      tol_ref = tol;

      if icase == 1
	nzk = 10;
      elseif icase == 2 || icase == 3
	nzk = 9;
      end
      
      for i=1:nzk
      	disp('Compute numerical solution')
	
      	tstart = tic;
      	osoutput = helmos(icase,tol,bctype,nsub,zk,npanvec,xylim,Ng,helmosopts);
      	dt = toc(tstart);
      	totaltime=[totaltime;dt];
      	unum = osoutput.unum;
      	np0=osoutput.np
      	np=[np;np0];
    
	disp(' ')
	disp('Compute reference solution')
	npanvec_ref = round(npanvec*1.5);
%	npanvec_ref = round(npanvec*2.0);
	osoutput_ref = helmos(icase,tol_ref,bctype,nsub_ref,zk, ...
			      npanvec_ref,xylim,Ng,helmosopts);
	uref = osoutput_ref.unum;
	
	zktot = [zktot;zk];
          
	rerr2 = norm(unum-uref)/norm(uref)
	relerr=[relerr;rerr2];
      
      %	rerrinf = norm(unum-uref,Inf)/norm(uref,Inf);
      %	rellinferr = [rellinferr;rerrinf];
          
	iter0 = osoutput.iter;
	iter = [iter; iter0];
          
	timeinfo = osoutput.timeinfo;
      
	buildtime0 = timeinfo(1)+timeinfo(2); % matrix correction + RCIP time
	buildtime = [buildtime;buildtime0];
      
	gmrestime0 = timeinfo(3);
	gmrestime = [gmrestime; gmrestime0];
%       double the wavenumber and npanvec
	zk = zk*2;
	npanvec = npanvec*2;
      end

      rellen = Ltot*zktot/(2*pi); % length relative to the wavelength
      [relerr rellen iter]
      [totaltime gmrestime buildtime np]

%      timing_plot(totaltime,gmrestime,buildtime,np,bctype,icase)
%      zk_iter_plot(iter,rellen,bctype,icase)
    end    

    
    if plotfield == 1
      disp('First solve: calculate the numerical solution')
      osoutput = helmos(icase,tol,bctype,nsub,zk,npanvec,xylim,Ng,helmosopts);
      unum = osoutput.unum;

      disp(' ')
      disp('Second solve: calculate the reference solution')
      npanvec_ref = round(npanvec*1.5);
      nsub_ref = round(nsub*1.5);
      osoutput_ref = helmos(icase,tol,bctype,nsub_ref,zk, ...
			    npanvec_ref,xylim,Ng,helmosopts);
      uref = osoutput_ref.unum;

      xg = osoutput.xg;
      yg = osoutput.yg;
      if icase >= 1
	fieldplot_abs(abs(unum),xg,yg,xylim,Ng,bctype,icase);
	field_error_plot(unum-uref,xg,yg,xylim,Ng,bctype,icase)
	disp(' ')
      else
	fieldplot_real(real(unum),xg,yg,xylim,Ng,bctype,icase);
	field_error_plot(real(unum-uref),xg,yg,xylim,Ng,bctype,icase)
	disp(' ')
      end
    end
  end  



  function osoutput = helmos(icase,tol,bctype,nsub,zk,npanvec,xylim,Ng,helmosopts)
% **************************************************************************
% driver for the Helmholtz open surface problem
%
% input parameters:
% icase: specifies the geometry 
%        icase=0: a line segment
%        icase=1: spiral
%        icase=2: corners
%        icase=3: triple junctions
%        icase=4:
% tol: desired precision for GMRES and FMM    
% bctype: = 'n': Neumann boundary condition
%         = 'd': Dirichlet boundary condition
% nsub: number of dyadic refinement toward the singular point on each edge
% zk: the wavenumber
% xylim: [xmin xmax ymin ymax] specifies the rectangle where the field is
%                              computed.
% Ng: the field is computed on NgxNg equispaced grid on the above rectangle
%
%
% output:
% osoutput - structure containing the following fields
% xg,yg: xy coordinates of the target points on the NgxNg tensor grid
% zsrc: source points on the boundary
% unum: numerical solution at NgxNg tensor grid
% iter: number of GMRES iterations
% residue: actual GMRES residue
% timeinfo: timing results for various steps in the whole algorithm
%    
%    
%    
%    
%    
%    
% **************************************************************************
% *** Solves the boundary value problem for the Helmholtz equation
% *** on piecewise smooth open curves.
% *** The boundary condition could be either Dirichlet or Neumann.
% ***
% *** A rectangular Cartesian grid of target points for field evaluation.
% *** 
% *** For the Dirichlet problem, the representation is u = S (-D') \sigma,
% *** where S is the Helmholtz single layer potential, and 
% *** D' is the normal derivative of the Helmholtz double layer potential.
% *** 
% *** For the Neumann problem, the representation is u = D(-S) \sigma,
% *** where D is the Helmholtz double layer potential.
% *** And the resulting boundary integral equation is D'(-S) \sigma = f.
% ***
% *** The integral equation formulations are from Bruno and Lintner.
% *** 
% *** The numerical method follows from the RCIP tutorial.
% ***
% *** modified from helmbvp5s_fast.m
% *** 
% *** Modifications: 1. added a parameter precond to specify different
% *** wavenumber for the preconditioner. Numerical experiments indicate
% *** that using the same wavenumber seems best. Somehow, D'_k S_{ik}
% *** doesn't work at all! Needs further investigation.
% *** 
% *** 2. changed the construction of sparse matrix for the quadrature
% *** correction part of the system matrix. This was NOT done correctly
% *** before and was very slow. The current version of osbuildmat_correction
% *** is at least 100 times faster than the one in helmbvp5s_fast.m, simply
% *** due to a better way of calling the function "sparse".
% *** 
% *** 3. changed the construction of sparse R matrices to efficiently handle
% *** the case of many identical singular points.
% *** 
% *** 
% *** Last modified on 2024-10-10
% **************************************************************************

   
    dlim=1.1;
    Ngtot=Ng*Ng;  % number of target points on the rectangular Cartesian grid
    xg=linspace(xylim(1),xylim(2),Ng);
    yg=linspace(xylim(3),xylim(4),Ng);
    zg=zeros(Ngtot,1);
    for k=1:Ng
      zg((k-1)*Ng+(1:Ng))=xg(k)+1i*yg;
    end

    timeinfo=[];
    ztarg=[];
    if isfield(helmosopts,'ztarg')
      ztarg = helmosopts.ztarg;
    end
%   obtain physical and geometric parameters
    osparams = ossetup(icase,helmosopts);
    src = [];
    if isfield(osparams,'src')
      src = osparams.src;
    end
    if isfield(osparams,'ncorner')
      ncorner = osparams.ncorner;
    end
    curvetype = [];
    if isfield(osparams,'curvetype')
      curvetype = osparams.curvetype;
    end
    if isfield(osparams,'cpars')
      cpars = osparams.cpars;
    end
    if isfield(osparams,'cparams')
      cparams = osparams.cparams;
    end
    if isfield(osparams,'ncurve')
      ncurve = osparams.ncurve;
    end
    if isfield(osparams,'nsing')
      nsing = osparams.nsing;
    end
    if isfield(osparams,'singinfo')
      singinfo = osparams.singinfo;
    end
    if isfield(osparams,'zsing')
      zsing = osparams.zsing;
    end

%   number of Gauss-Legendre nodes on each chunk
    ngl = 16;
    
%   this U is used with WfrakLinit that does matrix-vector
%   product instead of a linear solve.
%   U converts function values at Legendre nodes to the coefficients
%   of its Legendre expansion.
    [T16,W16,U,~]=legeexps(ngl);
%   logarithmic correction for self-interactions for each panel,
%   same for all panels. Hence the precomputation.    
    ML0=-log(abs(T16-T16'));
    ML0(1:ngl+1:end)=0;
    ML0=ML0+WfrakLinit(U,0,1,T16)./W16';
    
%   define functions for curves
    fcurve = cell(1,ncurve);
    for icurve=1:ncurve
      fcurve{icurve} = @(t) funcurve(t,icurve,cpars{icurve},icase,curvetype);
    end

%   npanvec
    np=sum(npanvec)*ngl;
    disp(['Total number of unknowns = ',num2str(sum(npanvec)*ngl)])

%   discretize the boundary
    chnkr = cell(1,ncurve);
    izind = cell(1,ncurve);   % indices on curves 
    ztot = [];
    for icurve = 1:ncurve
      [sinter,sinterdiff,~] = panelinit12(0,0,npanvec(icurve));
      
      chnkr{icurve} = chunkerfunc(fcurve{icurve},sinter,sinterdiff,T16,W16);
      if icase == 1 
        arclen=sum(chnkr{icurve}.awzp);
        npan=chnkr{icurve}.npan;
        [sinter,sinterdiff]=myremesh(arclen,T16,W16,npan,fcurve{icurve});
        chnkr{icurve}=chunkerfunc(fcurve{icurve},sinter,sinterdiff,T16,W16);
      end
      z = chnkr{icurve}.z;
      s = chnkr{icurve}.s;
      izind{icurve} = length(ztot)+(1:length(z));
      ztot = [ztot;z];  
    end

    chnkrtot=mergechnkr(chnkr);
    stot = chnkrtot.s;

% ************************************************************
% info needed for backward recursion to compute the density    
% on the fine mesh
% ************************************************************
    izindS = cell(1,nsing);   % indices on junctions
    for ising = 1:nsing
      tmp=[];
      clist  = singinfo{ising}.clist;
      nedge  = singinfo{ising}.nedge;	
      isstart= singinfo{ising}.isstart;		
      for i = 1:nedge
        if isstart(i)
          tmp=[tmp izind{clist(i)}(1:2*ngl)];
        else
          tmp=[tmp izind{clist(i)}(end-2*ngl+1:end)];
        end
      end
      izindS{ising}=tmp;
    end

    nsevec=zeros(1,nsing);   % number of stored levels 
    for ising = 1:nsing
      ds=min(abs(zg-zsing(ising)));

      clist=singinfo{ising}.clist;
      nedge=singinfo{ising}.nedge;
      isstart=singinfo{ising}.isstart;
      panlenvec=zeros(nedge,1);
      for icurve=1:nedge
        if isstart(icurve)
          panlenvec(icurve)=chnkr{clist(icurve)}.panlen(1);
        else
          panlenvec(icurve)=chnkr{clist(icurve)}.panlen(end);
        end
      end
      panlen=max(panlenvec);
        
%      premavec(ising)=nsub-ceil(log2(panlen/ds))-5;
      prema = nsub-ceil(log2(panlen/ds))-1;
      if prema >= nsub
        prema=nsub;
      end
      nsevec(ising)=nsub-prema+1;
    end

%   now deal with "identical" singular points
%   Use the largest nse values for all "identical" singular points    
    for ising = 1:nsing
      if singinfo{ising}.ifRcomp == 0
	isingcopy = singinfo{ising}.isingcopy;
	nse1 = nsevec(ising);
	nse2 = nsevec(isingcopy);
	if nse1>nse2
	  nsevec(isingcopy) = nse1;
	end
      end
    end

    nsevec
    
%     *** fine discretization ***
    chnkrfin = cell(1,ncurve);
    for icurve=1:ncurve
      nsub1=nsevec(cparams{icurve}.sing(1));
      nsub2=nsevec(cparams{icurve}.sing(2));
      sinter=chnkr{icurve}.sinter;
      sinterdiff=chnkr{icurve}.sinterdiff;
      [sinterfin,sinterdiffin,~]=panelinit12b(nsub1,nsub2,sinter,sinterdiff);
      chnkrfin{icurve}=chunkerfunc(fcurve{icurve},sinterfin, ...
				   sinterdiffin,T16,W16);
    end
      
%   build the system matrix
    rpars = [];
    rpars.ngl = ngl;
    rpars.zk = zk;
    rpars.T = T16;
    rpars.W = W16;
    rpars.U = U;
    rpars.ML0 = ML0;
    rpars.bctype=bctype;

%   0 -> preconditioner wavenumber = k; 1 -> preconditioner wavenumber = ck
    precond=0;
    if precond == 0
      zkS=zk;
      zkT=zk;
    end
%   experimenting preconditioners with ck as the wavenumber
%   numerical results show that (a) c=1i does NOT work at all
%   (b) c=1 seems best in terms of number of iterations in GMRES
%   (c) there is a small domain in the complex plane for c where
%   the preconditioner kind of works. But just as (a), there is 
%   a large domain where the preconditioner does not work
%   This is somewhat surprsing to me since one would think 
%   locally there is not much difference. 
    if precond == 1
      if bctype == 'd'
      	zkS = zk;
      	zkT = (1+0.5i)*zk;
      elseif bctype == 'n'
      	zkS = (1+0.5i)*zk;
      	zkT = zk;
      end
    end

    rpars.precond=precond; 
    rpars.zkS = zkS;
    rpars.zkT = zkT;

      
    disp(' ')
    disp('Step 1: build the sparse system matrix correction.')
    start = tic;
    [Sc,Tc] = osbuildmat_correction(chnkr,rpars);      
    dt = toc(start);
    disp(['System matrix construction time = ', num2str(dt), ' seconds'])
    timeinfo=[timeinfo;dt];


      
%   compute the preconditioner R in the RCIP method for all singular points
    disp(' ')
    disp('Step 2: compute the preconditioners for singular points.')
    start = tic;
  
    inds = [0; cumsum(npanvec)];

    R = cell(1,nsing);
    Rstor = cell(1,nsing);
    Kstor = cell(1,nsing);
    geomstor = cell(1,nsing);

    R1    = cell(1,nsing);
    R2    = cell(1,nsing);
    R3    = cell(1,nsing);
    R4    = cell(1,nsing);
    R13   = cell(1,nsing);
    R4213 = cell(1,nsing);
    Sstar = cell(1,nsing);
    Tstar = cell(1,nsing);

    rind=[];cind=[];rind1=(1:np)';
    for ising=1:nsing
      clist = singinfo{ising}.clist;
      isstart = singinfo{ising}.isstart;
      nedge = singinfo{ising}.nedge;

      if singinfo{ising}.ifRcomp == 1
        disp(['ising = ', num2str(ising)])
        rparslocal         = [];
        rparslocal.ngl     = ngl;
        rparslocal.zkS     = zkS;
        rparslocal.zkT     = zkT;
        rparslocal.T       = T16;
        rparslocal.W       = W16;
        rparslocal.U       = U;
        rparslocal.ML0     = ML0;
        rparslocal.bctype  = bctype;
        rparslocal.precond = precond;
	  
        cparslocal = cell(1,nedge);
        for i=1:nedge
          cparslocal{i} = cpars{clist(i)};
          cparslocal{i}.islocal = isstart(i);
        end

        fcurvelocal = cell(1,nedge);
        for i=1:nedge
          fcurvelocal{i} = @(t) funcurve(t,clist(i),cparslocal{i},icase,curvetype);
        end
%       panel size at the coarsest level in the parameter space
        h0 = zeros(1,nedge); % panel attached to the singular point
        h1 = zeros(1,nedge); % next panel
        for i=1:nedge
          if isstart(i)
            h0(i) = chnkr{clist(i)}.sinterdiff(1);
            h1(i) = chnkr{clist(i)}.sinterdiff(2);
          else
            h0(i) = chnkr{clist(i)}.sinterdiff(end);
            h1(i) = chnkr{clist(i)}.sinterdiff(end-1);
          end
        end
        h01 = [h0; h1]; % send both panel lengths to Rcomp
        
        ndim=2;
        rcipstruct = rcipsetup(ngl,ndim,nedge,isstart,T16,W16);
        geomstor{ising}.PbcL  = rcipstruct.PbcL;
        geomstor{ising}.starL = rcipstruct.starL;
        geomstor{ising}.circL = rcipstruct.circL;

      	[R{ising},Rstor{ising},Kstor{ising}] = Rcomp(nedge,ndim,nsub, ...
        		rcipstruct,h01,isstart,fcurvelocal, ...
      			nsevec(ising),rparslocal,T16,W16,ngl);

      	nR = length(R{ising})/ndim;
        R1{ising}    = R{ising}(1:nR,1:nR);
        R2{ising}    = R{ising}(nR+(1:nR),1:nR);
        R3{ising}    = R{ising}((1:nR),nR+(1:nR));
        R4{ising}    = R{ising}(nR+(1:nR),nR+(1:nR));
        R13{ising}   = R1{ising}\R3{ising};
        R4213{ising} = R4{ising}-R2{ising}*R13{ising};
      else
      	isingcopy = singinfo{ising}.isingcopy;
      	geomstor{ising} = geomstor{isingcopy};
      	R{ising} = R{isingcopy};
      	Rstor{ising} = Rstor{isingcopy};
      	Kstor{ising} = Kstor{isingcopy};
	
        R1{ising}    = R1{isingcopy};
        R2{ising}    = R2{isingcopy};
        R3{ising}    = R3{isingcopy};
        R4{ising}    = R4{isingcopy};
        R13{ising}   = R13{isingcopy};
        R4213{ising} = R4213{isingcopy};
      	Sstar{ising} = Sstar{isingcopy};
      	Tstar{ising} = Tstar{isingcopy};
      end
	
      starind = [];
      for i=1:nedge
        if isstart(i)
          starind = [starind inds(clist(i))*ngl+(1:2*ngl)];
        else
          starind = [starind inds(clist(i)+1)*ngl-fliplr(0:2*ngl-1)];
        end
      end

      if singinfo{ising}.ifRcomp == 1
      	[Sstar{ising},Tstar{ising}] = osbuildmat_smooth_starind(chnkrtot, ...
								starind,zkS,zkT);
      end
	
      starind=starind(:);
      rind1=setdiff(rind1,starind); % update nonzero diagonal entry indices
      ns=length(starind);
      rind0=repmat(starind,ns,1);
      cind0=repmat(starind',ns,1);cind0=cind0(:);

      rind=[rind;rind0];
      cind=[cind;cind0];
	
      Sc(starind,starind) = 0;
      Tc(starind,starind) = 0;
    end % end of for ising=1:nsing loop

    R1G0=[];R2G0=[];R13G0=[];R4213G0=[];Sss0=[];Tss0=[];
    for ising=1:nsing
      R1G0    = [R1G0;R1{ising}(:)];
      R2G0    = [R2G0;R2{ising}(:)];
      R13G0   = [R13G0;R13{ising}(:)];
      R4213G0 = [R4213G0;R4213{ising}(:)];
      Sss0    = [Sss0;Sstar{ising}(:)];
      Tss0    = [Tss0;Tstar{ising}(:)];
    end

    v = ones(length(rind1),1);
    IG = sparse(rind1,rind1,v,np,np);
    
    R1G0 = [R1G0;v]; R4213G0 = [R4213G0;v];
    rind2 = [rind;rind1]; cind2 = [cind;rind1];
    
    R1G    = sparse(rind2,cind2,R1G0,np,np);
    R4213G = sparse(rind2,cind2,R4213G0,np,np);

    R2G    = sparse(rind,cind,R2G0,np,np);
    R13G   = sparse(rind,cind,R13G0,np,np);

    Sss    = sparse(rind,cind,Sss0,np,np);
    Tss    = sparse(rind,cind,Tss0,np,np);

    dt = toc(start);
    disp(['Preconditioner construction time = ', num2str(dt), ' seconds'])
    timeinfo=[timeinfo;dt];



      
%   right hand side: irhs = 0 -> monopole; irhs=1 -> planewave
    irhs=1;
    if ~isempty(ztarg)
      irhs=2;
    end

    if icase == 4;
%      irhs = 0;
    end
    
    if irhs == 0 % monopole
      if bctype == 'n' % Neumann
        zdiff=ztot-src;
        tmp=zk*abs(zdiff);
        nz=chnkrtot.nz;
        rhs=tmp.*besselh(1,tmp).*real(nz./zdiff);
      else % Dirichlet
        rhs=-besselh(0,1,zk*abs(ztot-src)); 
      end
    elseif irhs == 1 % plane wave
      if icase == 1
	if bctype == 'd'
	  theta = 3*pi/4; % incident angle
	end
	if bctype == 'n'
	  theta = pi; % incident angle
	end
      elseif icase == 2
	if bctype == 'd'
	  theta =  pi/4;
	end
	if bctype == 'n'
	  theta = -pi/4;
	end
      elseif icase == 3
	if bctype == 'd'
	  theta = -2*pi/3;
	end
	if bctype == 'n'
	  theta = pi/3;
	end
      elseif icase == 4
	if bctype == 'd'
	  theta = pi/6;
	end
	if bctype == 'n'
	  theta = 4*pi/3;
	end
      end
      dd=exp(1i*theta);

      if bctype == 'd' % Dirichlet
	rhs=-exp(1i*abs(zk)*real(ztot*conj(dd)));
      end
      
      if bctype == 'n' % Neumann
	rhs=-exp(1i*abs(zk)*real(ztot*conj(dd)));
        rhs=1i*zk*rhs.*real(chnkrtot.nz*conj(dd));
      end
    elseif irhs == 2
      x=stot*2-1;
      rhs=4*x.^3+2*x.^2-3*x-1;
    end



      
%   solve the linear system using gmres
    disp(' ')
    disp('Step 3: solve the linear system via GMRES.')
    start = tic;
      
%   use reduced equation, similar to (60) in RCIP tutorial
    [rhotildetot,iter] = myGMRESR2(zkS,zkT,chnkrtot,Sss,Tss,Sc,Tc,IG,R1G,R2G, ...
				   R13G,R4213G,rhs,np,3000,tol,bctype);
    
%    [rhotildetot,iter] = myGMRESR3(zkS,zkT,chnkrtot,Sss,Tss,Sc,Tc,IG,R1G,R2G, ...
%				   R13G,R4213G,rhs,np,150,tol,bctype,10000);

    
%    max_iter = 2000;
%    restart = 100;
%    x0 = zeros(np,1);
%    keep = restart/2;
%    [rhotildetot,flag,relres,iter] = gmres_restart_augmented(zkS,zkT,chnkrtot, ...
%				Sss,Tss,Sc,Tc,IG,R1G,R2G,R13G,R4213G,bctype, ...
%				rhs,tol,max_iter,restart,x0,keep);

    disp(['GMRES iterations = ',num2str(iter)])
    rho1hattot=R1G*rhotildetot;

    if bctype == 'd' % B=-T
      rho2tildetot = Tss*rho1hattot-Toper_fmm(tol,zkT,chnkrtot,rho1hattot)-Tc*rho1hattot;
    elseif bctype == 'n' % B=-S
      rho2tildetot = Sss*rho1hattot-Soper_fmm(tol,zkS,chnkrtot,rho1hattot)-Sc*rho1hattot;
    end
    rho2hattot = R2G*rhotildetot+R4213G*rho2tildetot;
    rho1tildetot = rhotildetot-R13G*rho2tildetot;
    dt = toc(start);
    disp(['Time on GMRES = ', num2str(dt), ' seconds'])
    disp(' ')
    timeinfo=[timeinfo;dt];




    
    disp(' ')
    disp('Step 4: reconstruct density on the fine mesh.')
% *** rho2hat on Curves ***      
    tstart=tic;
    rho2hatC=cell(1,ncurve);
    for icurve = 1:ncurve
      rho2hatC{icurve} = rho2hattot(izind{icurve});
    end
 
% *** reconstructed refined rho2hat in Singular zones ***
    rho2hatSfin=myreconstruct(Kstor,Rstor,geomstor,rho1tildetot, ...
			      rho2tildetot,izindS,nsevec,singinfo,nsing,ngl,ndim);
      
% *** reconstructed refined rho2hat on Curves ***      
    rho2hatCfin=SfintoCfin(rho2hatSfin,rho2hatC,cparams,singinfo, ...
			   nsevec,npanvec,ncurve,ngl);
    rho2hatfintot=[];
    for icurve=1:ncurve
      rho2hatfintot=[rho2hatfintot; rho2hatCfin{icurve}];
    end

    chnkrfintot=mergechnkr(chnkrfin);
    dt = toc(tstart);
    disp(['Fine density construction time = ', num2str(dt), ' seconds'])
    timeinfo=[timeinfo;dt];

    disp(' ')
    disp('Step 5: evaluate the field.')
    tstart=tic;
    if ~isempty(ztarg)
      unum=Starg(ztarg,ztot,zk)*(rho2hattot.*chnkrtot.awzp);
    else
      [unum,outlist]=rectgrid_fieldcomp_fmm(tol,zk,chnkrfintot, ...
			rho2hatfintot,xg,yg,zg,dlim,cpars,ngl,icase,bctype);
% compute the total field if the incident field is a plane wave
      if irhs == 1
      	unum = unum+sum(exp(1i*zk*real(zg*conj(dd))),2);
      elseif irhs == 0
	unum = unum+besselh(0,1,zk*abs(zg-src)); 
      end
    end
    dt = toc(tstart);
    disp(['Time on fieldcomp = ', num2str(dt), ' seconds'])
    timeinfo=[timeinfo;dt];

    osoutput=[];
    osoutput.xg       = xg;
    osoutput.yg       = yg;
    osoutput.zsrc     = ztot;
    osoutput.unum     = unum;
    osoutput.iter     = iter;
    osoutput.np       = np;
    osoutput.timeinfo = timeinfo;
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Backward recursion subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function rho2hatSfin=myreconstruct(Kstor,Rstor,geomstor,rho1tildetot, ...
		rho2tildetot,izindS,nsevec,singinfo,nsing,ngl,ndim)
% *** Reconstruction of rhohatfin in Singular zones ***
    rhohatSfin=cell(1,nsing);   % all densities
    rho2hatSfin=cell(1,nsing);  % density of interest    
    for ising = 1:nsing
      nedge=singinfo{ising}.nedge;
      isstart=singinfo{ising}.isstart;
      nse=nsevec(ising);
      pts2=ngl*(nse+2);   % number of points on an edge
      rhohatSfin{ising}=zeros(ndim*nedge*pts2,1);      
      rhot=[rho1tildetot(izindS{ising});rho2tildetot(izindS{ising})];	
      PbcL =geomstor{ising}.PbcL;
      starL=geomstor{ising}.starL;
      circL=geomstor{ising}.circL;
      nsys=3*ngl*nedge*ndim;
      nR=2*ngl*nedge*ndim;
      for k = nse:-1:1
        Kcirc=Kstor{ising}(:,:,nse-k+1);
        Kcirc(starL,starL)=zeros(nR);
        MAT=eye(nsys)+Kcirc;
        tmp=PbcL*rhot;
        rhot=tmp-Kcirc*SchurBanaLs(MAT,Rstor{ising}(:,:,nse-k+1),tmp,starL,circL);
        myind1=[];
        for i = 1:nedge
          if isstart(i)
            myindL=(k+1)*ngl+(1:ngl);
            myind1=[myind1 myindL+(i-1)*pts2];
          else
            myindL=(nse-k)*ngl+(1:ngl);
            myind1=[myind1 myindL+(i-1)*pts2];
          end
        end
        myind=myind1;
        for i = 1:ndim-1
          myind=[myind myind1+i*nedge*pts2];
        end
        rhohatSfin{ising}(myind)=rhot(circL);
        rhot=rhot(starL);
      end
      myind1=[];
      for i = 1:nedge
        if isstart(i)
          myindL=1:2*ngl;
          myind1=[myind1 myindL+(i-1)*pts2];
        else
          myindL=pts2-2*ngl+1:pts2;
          myind1=[myind1 myindL+(i-1)*pts2];
        end
      end
      myind=myind1;
      for i = 1:ndim-1
        myind=[myind myind1+i*nedge*pts2];
      end
      rhohatSfin{ising}(myind)=Rstor{ising}(:,:,nse)*rhot;
%     *** we are only interested in the second density ***
      rho2hatSfin{ising}=rhohatSfin{ising}(nedge*pts2+1:2*nedge*pts2);
    end      
  end

 
  
  function rho2hatCfin=SfintoCfin(rho2hatSfin,rho2hatC,cparams,singinfo, ...
				  nsevec,npanvec,ncurve,ngl)
    rho2hatCfin=cell(1,ncurve);
    for icurve = 1:ncurve
      ising1=cparams{icurve}.sing(1);
      ising2=cparams{icurve}.sing(2);
      nse1=nsevec(ising1);
      nse2=nsevec(ising2);
      clist1=singinfo{ising1}.clist;
      clist2=singinfo{ising2}.clist;
      iedge1=find(clist1==icurve);
      iedge2=find(clist2==icurve);
      pts21=ngl*(nse1+2);
      pts22=ngl*(nse2+2);
      rho2hatCfin{icurve}= ...
	    [rho2hatSfin{ising1}((iedge1-1)*pts21+1:iedge1*pts21);
	     rho2hatC{icurve}(2*ngl+1:(npanvec(icurve)-2)*ngl);
	     rho2hatSfin{ising2}((iedge2-1)*pts22+1:iedge2*pts22)];
    end
  end
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Field comp subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  function [unum,outlist]=fieldcomp_fmm(tol,omega,chnkr, ...
    rho2hat,zg,dlim,cpars,ngl,icase,bctype)
% vectorized version
    disp('Screening starts')
    start=tic;
%   *** mylist1 => some close evaluation on fine grid
    Ngtot=length(zg);

    npan   = chnkr.npan;
    zsrc      = chnkr.z;
    nz     = chnkr.nz;
    panlen = chnkr.panlen;
    zinter = chnkr.zinter;
    wzp    = chnkr.wzp;
    awzp   = chnkr.awzp;
    
    outlist=ones(Ngtot,1);

    mylist=outlist;
    mylist1=zeros(Ngtot,1);
    panlen2=panlen.^2;
    dlim2=dlim^2;
    targpan=cell(1,npan);
    nspec=0;
    for k=1:Ngtot
      if outlist(k)
        zdiff=zg(k)-zsrc;
        xd=real(zdiff);
        yd=imag(zdiff);
        d2=xd.*xd+yd.*yd; % sqrt is expensive to compute, so avoid it.
        d2=reshape(d2,ngl,npan);
        d=min(d2)'./panlen2;
        tmp=d<dlim2;
        if nnz(tmp)>0
      	  mylist1(k)=1;
          panno=find(tmp)';
          for kk=panno
            targpan{kk}=[targpan{kk} k];
          end
      	  nspec=nspec+1;
        end
      end
    end
    disp(['Number of evaluation points   = ',num2str(Ngtot)])
    dt = toc(start);
    disp(['Screening and sorting time    = ', num2str(dt), ' seconds'])
%    
    unum=zeros(Ngtot,1);
    start=tic;
%   use the FMM for regular evaluation    
    if bctype == 'd'
      unum=Starg_fmm(tol,omega,zg,zsrc,awzp,rho2hat);      
    elseif bctype == 'n'
      unum=Ktarg_fmm(tol,omega,zg,zsrc,nz,awzp,rho2hat);
    end
    dt = toc(start);
    disp(['Regular evaluation time       = ', num2str(dt), ' seconds'])
%
    disp(['Number of close-evaluation points =',num2str(nspec)])
    start=tic;
    sol=awzp.*rho2hat;
    for kk=1:npan
      k=targpan{kk};
      if ~isempty(k)
        a=zinter(1,kk);
        b=zinter(2,kk);
        myind=(kk-1)*ngl+1:kk*ngl;
	
        nz0=nz(myind);
        zsrc0=zsrc(myind);
        awzp0=awzp(myind);
        wzp0=wzp(myind);
        sol0=sol(myind);

      	ztrg=zg(k);
        U=wLCHS_src_precomp(a,b,zsrc0);

      	if bctype == 'n'
          [LogC,CauC]=wLCHS_target(a,b,ztrg,zsrc0,nz0,wzp0,awzp0,U,0);
          M=KtargCloseCorrection(ztrg,zsrc0,omega,wzp0,LogC,CauC);
          unum(k)=unum(k)+M*rho2hat(myind);
        else
          LogC=wLCHS_target(a,b,ztrg,zsrc0,nz0,wzp0,awzp0,U,0);
          M=StargCloseCorrection(ztrg,zsrc0,omega,LogC);
          unum(k)=unum(k)+M*sol0;
        end
      end
    end
    dt = toc(start);
    disp(['List1 correction time         = ', num2str(dt), ' seconds'])
  end


  
  function [unum,outlist]=rectgrid_fieldcomp_fmm(tol,omega,chnkr, ...
	rho2hat,xg,yg,zg,dlim,cpars,ngl,icase,bctype)
% field evaluation on rectangular grid    
    disp('Screening starts')
    start=tic;
%   *** mylist1 => some close evaluation on fine grid
    Ngtot=length(zg);

    npan   = chnkr.npan;
    zsrc   = chnkr.z;
    nz     = chnkr.nz;
    panlen = chnkr.panlen;
    zinter = chnkr.zinter;
    wzp    = chnkr.wzp;
    awzp   = chnkr.awzp;
    
    outlist=ones(Ngtot,1);

    targpan=getlist1_rectangulargrid(chnkr,xg,yg,dlim);
    nspec=0;
    for kk=1:npan
      nspec=nspec+length(targpan{kk});
    end
    
    disp(['Number of evaluation points   = ', num2str(Ngtot)])
    dt = toc(start);
    disp(['Screening and sorting time    = ', num2str(dt), ' seconds'])
%    
    unum=zeros(Ngtot,1);
    start=tic;
%   use the FMM for regular evaluation    
    if bctype=='d'
      unum=Starg_fmm(tol,omega,zg,zsrc,awzp,rho2hat);      
    elseif bctype=='n'
      unum=Ktarg_fmm(tol,omega,zg,zsrc,nz,awzp,rho2hat);
    else

    end
    dt = toc(start);
    disp(['Regular evaluation time       = ', num2str(dt), ' seconds'])
%
    disp(['Number of close-evaluation points = ',num2str(nspec)])
    start=tic;
    sol=awzp.*rho2hat;
    for kk=1:npan
      k=targpan{kk};
      if ~isempty(k)
        a=zinter(1,kk);
        b=zinter(2,kk);
        myind=(kk-1)*ngl+1:kk*ngl;
	
        nz0=nz(myind);
        zsrc0=zsrc(myind);
        awzp0=awzp(myind);
        wzp0=wzp(myind);
        sol0=sol(myind);

      	ztrg=zg(k);
        U=wLCHS_src_precomp(a,b,zsrc0);

      	if bctype=='n'
          [LogC,CauC]=wLCHS_target(a,b,ztrg,zsrc0,nz0,wzp0,awzp0,U,0);
          M=KtargCloseCorrection(ztrg,zsrc0,omega,wzp0,LogC,CauC);
          unum(k)=unum(k)+M*rho2hat(myind);
        else
          LogC=wLCHS_target(a,b,ztrg,zsrc0,nz0,wzp0,awzp0,U,0);
          M=StargCloseCorrection(ztrg,zsrc0,omega,LogC);
          unum(k)=unum(k)+M*sol0;
        end
      end
    end
    dt = toc(start);
    disp(['List1 correction time         = ', num2str(dt), ' seconds'])
  end




  function targpan=getlist1_rectangulargrid(chnkr,xg,yg,dlim)
% return list of targets for each source panel in chnkr     
% using simple bin sort on rectangular grid targets

    xmin=min(xg);
    xmax=max(xg);
    ymin=min(yg);
    ymax=max(yg);
    Nx=length(xg);
    Ny=length(yg);
    hx=(xmax-xmin)/(Nx-1);
    hy=(ymax-ymin)/(Ny-1);
    
    npan   = chnkr.npan;
    zsrc   = chnkr.z;
    panlen = chnkr.panlen;
    ngl    = chnkr.ngl;

    dlim2=dlim^2;
    
    targpan=cell(1,npan);

    for kk=1:npan
      myind=(kk-1)*ngl+1:kk*ngl;
      z=zsrc(myind);
      dd=dlim*panlen(kk);
      x1=min(real(z))-dd;
      x2=max(real(z))+dd;
      y1=min(imag(z))-dd;
      y2=max(imag(z))+dd;

      indx1=floor((x1-xmin)/hx)+1;
      indx1=max(indx1,1);
      indx2=ceil((x2-xmin)/hx)+1;
      indx2=min(indx2,Nx);
      nx=indx2-indx1+1;

      indy1=floor((y1-ymin)/hy)+1;
      indy1=max(indy1,1);
      indy2=ceil((y2-ymin)/hy)+1;
      indy2=min(indy2,Ny);
      ny=indy2-indy1+1;

      [xt,yt]=meshgrid(xg(indx1:indx2),yg(indy1:indy2));
      zt=xt+1i*yt;
      zt=zt(:);

      zdiff=zt-z.';
      xd=real(zdiff);
      yd=imag(zdiff);
      d2=xd.*xd+yd.*yd;
      d=min(d2,[],2)/panlen(kk)^2;
      d=reshape(d,ny,nx);
      [iy,ix]=find(d<dlim2);

      targpan{kk}= (iy'+indy1-1)+((ix'+indx1-1)-1)*Ny;
    end
  end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Operator evaluation subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function u=Starg(ztrg,zsrc,zk)
% single layer potential kernel
    u=0.5i*besselh(0,zk*abs(ztrg-zsrc.'));
  end


  
  function u=Starg_fmm(eps,zk,ztrg,zsrc,awzp,density)
% 
    sources=[real(zsrc).';imag(zsrc).'];
    srcinfo=[];
    srcinfo.sources=sources;
    srcinfo.charges=(awzp.*density).';

    targ=[real(ztrg).';imag(ztrg).'];

    pg=0;
    pgt=1;
    U=hfmm2d(eps,zk,srcinfo,pg,targ,pgt);
    u=2*U.pottarg;
    u=u(:);
  end

  
  
  function u=Soper_fmm(eps,zk,chnkr,density)
% single layer potential kernel
    z=chnkr.z;
    awzp=chnkr.awzp;
    srcinfo=[];
    srcinfo.sources=[real(z).';imag(z).'];
    srcinfo.charges=(awzp.*density).';
    
    pg=1;
    U1=hfmm2d(eps,zk,srcinfo,pg);
    u=2*U1.pot;
    u=u(:);
  end


  
  function u=Toper_fmm(eps,zk,chnkr,density)
% D'
    z=chnkr.z;
    nz=chnkr.nz;
    awzp=chnkr.awzp;
    sources=[real(z).';imag(z).'];    
    dipvec=[real(nz).';imag(nz).'];
    dipstr = density.*awzp;
    dipstr = dipstr(:).';
    srcinfo=[];
    srcinfo.sources=sources;
    srcinfo.dipvec=dipvec;
    srcinfo.dipstr=dipstr;

    targ=sources;
    pg=0;
    pgt=2;
    U1=hfmm2d(eps,zk,srcinfo,pg,targ,pgt);
    u=2*sum(dipvec.*U1.gradtarg,1);
    u=u(:);
  end

  
  
  function u=StargCloseCorrection(ztrg,zsrc,zk,LogC)
% single layer potential close evaluation correction    
    u=-besselj(0,zk*abs(ztrg-zsrc.')).*LogC/pi; 
  end
  
  
  
  function u=Ktarg(ztrg,zsrc,zk,wzp)
% double layer potential kernel times smooth Gauss-Legendre weights
    zdiff=ztrg-zsrc.';
    tmp=zk*abs(zdiff);
    u=0.5i*tmp.*besselh(1,tmp).*imag(wzp.'./zdiff);
  end



  function u=Ktarg_fmm(eps,zk,ztrg,zsrc,nz,awzp,density)
% double layer potential kernel times smooth Gauss-Legendre weights
    sources=[real(zsrc).';imag(zsrc).'];
    dipvec=[real(nz).';imag(nz).'];
    dipstr = density.*awzp;
    dipstr = dipstr(:).';
    srcinfo=[];
    srcinfo.sources=sources;
    srcinfo.dipvec=dipvec;
    srcinfo.dipstr=dipstr;

    targ=[real(ztrg).';imag(ztrg).'];

    pg=0;
    pgt=1;
    U=hfmm2d(eps,zk,srcinfo,pg,targ,pgt);
    u=2*U.pottarg;
    u=u(:);
  end

  
  
  function u=KtargCloseCorrection(ztrg,zsrc,zk,wzp,LogC,CauC)
% double layer potential close evaluation correction
    zdiff=ztrg-zsrc.';
    tmp=zk*abs(zdiff);
    u=tmp.*besselj(1,tmp).*imag(wzp.'./zdiff);
    u=-1/pi*(u.*LogC+real(CauC));   
  end
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  almost equispaced remesh in the arclength parameter
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [sinter,sinterdiff]=myremesh(arclen,T16,W16,npan,fcurve)
    panlen=arclen/npan;
    spanlen=1/npan;
    sinter=zeros(npan+1,1);
    for k=1:npan
      x0=sinter(k);
      aa=x0+0.1*spanlen;
      bb=min(x0+10*spanlen,1);
      for kk=1:55
        cc=(aa+bb)/2;
        awzpA=zinitA([x0 cc],T16,W16,fcurve);
        fA=sum(awzpA)-panlen;
        if fA>0
          bb=cc;
        else
          aa=cc;
        end
      end
      sinter(k+1)=cc;
    end
    sinter(2)=(sinter(1)+sinter(3))/2;
    sinter(npan)=(sinter(npan-1)+sinter(npan+1))/2;
    sinterdiff=diff(sinter);
  end


  
  function awzp=zinitA(sinter,T,W,fcurve)
    sdif=(sinter(2)-sinter(1))/2;
    s=(sinter(1)+sinter(2))/2+sdif*T;
    w=W*sdif;
    [~,zp]=fcurve(s);
    awzp=w.*abs(zp);
  end


  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Geometry setup subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function osparams = ossetup(icase,helmosopts)
    if icase == 0
      osparams = ossetup_linesegment();
    elseif icase == 1
      osparams = ossetup_spiral();
    elseif icase == 2
      ncorner = helmosopts.ncorner;
      osparams = ossetup_corners(ncorner);
    elseif icase == 3
      max_depth = helmosopts.max_depth;
      osparams = ossetup_branches(max_depth);
    elseif icase == 4
      osparams = ossetup_maze();
    end
  end

  
  
  function osparams = ossetup_linesegment()
    osparams = [];
%   number of curves
    ncurve = 1; 
%   number of singular points
    nsing=2;
    
%   endpoints   
    z=zeros(nsing,1); 
    z(1) =-1-0.2i;
    z(2) = 1-0.2i;

%   specify the boundary
    cparams = cell(1,ncurve);
    for icurve = 1:ncurve
      cparams{icurve}.ta = -1; % starting parameter value
      cparams{icurve}.tb = 1; % end parameter value
      cparams{icurve}.ifclosed = false; % whether each curve is closed or not
    end
    cparams{1}.sing=[1 2];
    
    cpars = cell(1,ncurve);
%   starting and end points of each line segment  
    for icurve = 1:ncurve
      cpars{icurve}.v0=z(cparams{icurve}.sing(1));
      cpars{icurve}.v1=z(cparams{icurve}.sing(2));
    end

%   structure contain information needed for singular points  
    singinfo = cell(1,nsing);
    for ising = 1:nsing
      singinfo{ising} = [];
      ifRcomp=1;
      isingcopy=0;
      if ising==1 % left endpoint
%       curve ID starting or ending at the singular point
        clist = [1];
%       indicators on whether each curve starts or ends at the singular point      
%       1 -- singular point is the starting point; 0 -- singular point is the end point   
        isstart = [1];
%       number of edges      
        nedge = 1;
      elseif ising == 2 % right endpoint
        clist = [1];
        isstart = [0];
        nedge = 1;
      end
  
      singinfo{ising}.clist = clist;
      singinfo{ising}.isstart = isstart;
      singinfo{ising}.nedge = nedge;
      singinfo{ising}.ifRcomp = ifRcomp;
      singinfo{ising}.isingcopy = isingcopy;
    end
  
    osparams.cpars = cpars;
    osparams.cparams = cparams;
    osparams.ncurve = ncurve;
    osparams.nsing = nsing;
    osparams.singinfo = singinfo;
    osparams.zsing = z;
  end


  
  function osparams = ossetup_spiral()
    osparams = [];
%   number of curves
    ncurve = 1; 
%   number of singular points
    nsing=2;
    
%   endpoints   
    z=zeros(nsing,1); 
    gammas=exp(-1-5i);
    gammae=exp( 1+5i);
    z(1)=gammas;
    z(2)=gammae;

%   specify the boundary
    cparams = cell(1,ncurve);
    for icurve = 1:ncurve
      cparams{icurve}.ta = 0; % starting parameter value
      cparams{icurve}.tb = 1; % end parameter value
      cparams{icurve}.ifclosed = false; % whether each curve is closed or not
    end
    cparams{1}.sing=[1 2];
    
    cpars = cell(1,ncurve);
%   starting and end points of each line segment  
    for icurve = 1:ncurve
      cpars{icurve}.v0=z(cparams{icurve}.sing(1));
      cpars{icurve}.v1=z(cparams{icurve}.sing(2));

      cpars{icurve}.ta=cparams{icurve}.ta;
      cpars{icurve}.tb=cparams{icurve}.tb;
      cpars{icurve}.a=exp(-(1+5i));
      cpars{icurve}.b=-1i*(2+10i);
    end

%   structure contain information needed for singular points  
    singinfo = cell(1,nsing);
    for ising = 1:nsing
      singinfo{ising} = [];
      ifRcomp=1;
      isingcopy=0;
      if ising==1 % left endpoint
        clist = [1];
        isstart = [1];
        nedge = 1;
      elseif ising == 2 % right endpoint
        clist = [1];
        isstart = [0];
        nedge = 1;
      end
  
      singinfo{ising}.clist = clist;
      singinfo{ising}.isstart = isstart;
      singinfo{ising}.nedge = nedge;
      singinfo{ising}.ifRcomp = ifRcomp;
      singinfo{ising}.isingcopy = isingcopy;
    end

    osparams.cpars = cpars;
    osparams.cparams = cparams;
    osparams.ncurve = ncurve;
    osparams.nsing = nsing;
    osparams.singinfo = singinfo;
    osparams.zsing = z;
  end


  

  function osparams = ossetup_corners(ncorner,curvetype)
% input parameters:
%    curvetype == 1 -> line segments
%    curvetype == 2 -> parabola + cosine curve
%    ncorner -> number of corners


    osparams = [];

    if nargin < 2
      curvetype = 2;
    end
    
%   number of corners
    if nargin < 1
      ncorner = 1;
    end
%   number of curves
    ncurve = ncorner+1; 
%   number of singular points
    nsing = ncorner+2;
    
    if curvetype == 1
      d=1;
%     endpoints   
      z=zeros(nsing,1); 
      z(1:2:ncorner+2)=(-ncorner/2-1)*d-0.5i+(1:2:ncorner+2)*d;
      z(2:2:ncorner+2)=(-ncorner/2-1)*d+(2:2:ncorner+2)*d+0.5i;
    elseif curvetype == 2
      d=(1+pi/2)/2;
      z=zeros(nsing,1);
      z(1:2:ncorner+2) = ((1:2:(ncorner+2))-1)*d;
      z(2:2:ncorner+2) = ((2:2:(ncorner+2))-2)*d+1+1i;
    end

    if curvetype == 2
      a=1;b=0;c1=1;c2=0;x0=0; % parabola
      z0=1;z1=1; % cosine curve
    end
    
%   specify the boundary
    cparams = cell(1,ncurve);
    for icurve = 1:ncurve
      cparams{icurve}.ta = 0; % starting parameter value
      cparams{icurve}.tb = 1; % end parameter value
      cparams{icurve}.ifclosed = false; % whether each curve is closed or not
    end

    for icurve = 1:ncurve
      cparams{icurve}.sing = [icurve icurve+1];
    end
    
    cpars = cell(1,ncurve);
%   starting and end points of each line segment
    if curvetype == 1
      for icurve = 1:ncurve
	cpars{icurve}.v0=z(cparams{icurve}.sing(1));
	cpars{icurve}.v1=z(cparams{icurve}.sing(2));    
      end
    elseif curvetype == 2
      for icurve = 1:2:ncurve
	cpars{icurve}.a  = a;
	cpars{icurve}.b  = b;
	cpars{icurve}.x0 = real(z(icurve));
	cpars{icurve}.c1 = c1;
	cpars{icurve}.c2 = real(z(icurve));
      end
      for icurve = 2:2:ncurve
	cpars{icurve}.z0 = real(z(icurve));
	cpars{icurve}.z1 = z1;
      end
    end
%   structure contain information needed for singular points  
    singinfo = cell(1,nsing);
    for ising = 1:nsing
      singinfo{ising} = [];
      ifRcomp=1;
      isingcopy=0;
      if ising==1 % left endpoint
        clist = [1];
        isstart = [1];
        nedge = 1;
      elseif ising == ncorner+2 % right endpoint
        clist = [ncurve];
        isstart = [0];
        nedge = 1;
      else % corners
        clist = [ising-1,ising];
        isstart = [0,1];
        nedge = 2;
      end

      if ising >= 4 && ising < ncorner+2
	ifRcomp=0;
	isingcopy=mod(ising,2)+2;
      end
      
      singinfo{ising}.clist = clist;
      singinfo{ising}.isstart = isstart;
      singinfo{ising}.nedge = nedge;
      singinfo{ising}.ifRcomp = ifRcomp;
      singinfo{ising}.isingcopy = isingcopy;
    end

    osparams.cpars = cpars;
    osparams.cparams = cparams;
    osparams.ncurve = ncurve;
    osparams.curvetype = curvetype;
    osparams.ncorner = ncorner;
    osparams.nsing = nsing;
    osparams.singinfo = singinfo;
    osparams.zsing = z;
  end

  
  function osparams = ossetup_branches(max_depth)
% input parameters:
%    max_depth -> level of recursive construction. 0 -> one branch point

    osparams = [];

    if nargin < 1
      max_depth = 0;
    end
%   sc -> scaling factor for the branch_length on the next level
    sc = 0.8;
    
%   number of branch points
    nbranch = 2^(max_depth+1)-1

%   number of curves
    ncurve = 2^(max_depth+2)-1;
    
%   number of singular points
    nsing = ncurve+1;
    
%   endpoints   
%   specify the boundary
    cparams = cell(1,ncurve);
    for icurve = 1:ncurve
      cparams{icurve}.ta = 0; % starting parameter value
      cparams{icurve}.tb = 1; % end parameter value
      cparams{icurve}.ifclosed = false; % whether each curve is closed or not
    end

    angle = zeros(max_depth+1);
    angle(1) = pi/4; 
    for i=1:max_depth
      angle(i+1)=angle(i)/1.1; % change 1.1 to 1 to make all branch points the same
    end

    z = zeros(nsing,1);z(1)=0;
    dz = zeros(nsing-1,1);dz(1)=1i*10;
    z(2)=z(1)+dz(1);

    ind = []; ind=[ind; 1];ie=ind(end);
    icurve=1;
    cparams{icurve}.sing = [icurve icurve+1];    
    for i=0:max_depth
      ns = length(ind);
      ie = ind(ns);
      dz(ie+(1:ns)*2-1) = dz(ind)*sc*exp(-1i*angle(i+1));
      dz(ie+(1:ns)*2) = dz(ind)*sc*exp(1i*angle(i+1));
      z(ie+(1:ns)*2) = z(ind+1)+dz(ie+(1:ns)*2-1);
      z(ie+(1:ns)*2+1) = z(ind+1)+dz(ie+(1:ns)*2);
      for k=1:ns
	icurve=icurve+1;
	cparams{icurve}.sing = [ind(k)+1, ie+2*k];    
	icurve=icurve+1;
	cparams{icurve}.sing = [ind(k)+1, ie+2*k+1];    
      end

      ind = ie+(1:2*ns);
    end

    cpars = cell(1,ncurve);
%   starting and end points of each line segment
    for icurve = 1:ncurve
      cpars{icurve}.v0=z(cparams{icurve}.sing(1));
      cpars{icurve}.v1=z(cparams{icurve}.sing(2));
    end
    
%   structure contain information needed for singular points  
    singinfo = cell(1,nsing);
    for ising = 1:nsing
      singinfo{ising} = [];
      ifRcomp=1;
      isingcopy=0;
      if ising==1 % root endpoint
        clist = [1];
        isstart = [1];
        nedge = 1;
      elseif ising>=nsing-2^(max_depth+1)+1 % leaf endpoints
        clist = [ising-1];
        isstart = [0];
        nedge = 1;
      else % branch points
	c1 = ising-1;
        clist = [c1,2*c1,2*c1+1];
        isstart = [0,1,1];
        nedge = 3;
      end

      if ising>nsing-2^(max_depth+1)+1 % leaf endpoints
	ifRcomp = 0;
	isingcopy = nsing-2^(max_depth+1)+1;
      end
      singinfo{ising}.clist = clist;
      singinfo{ising}.isstart = isstart;
      singinfo{ising}.nedge = nedge;
      singinfo{ising}.ifRcomp = ifRcomp;
      singinfo{ising}.isingcopy = isingcopy;
    end

    osparams.cpars = cpars;
    osparams.cparams = cparams;
    osparams.ncurve = ncurve;
    osparams.nsing = nsing;
    osparams.singinfo = singinfo;
    osparams.zsing = z;
  end

  
 
  function osparams = ossetup_maze(M,N)
% input parameters:
%    draw an MxN maze

    osparams = [];

    if nargin < 1
      M=10;
      N=10;
    end

    persistent maze ifplot src
    
    if isempty(maze)
%      maze = generateMaze(M,N);
      maze0 = load("maze1.data");
      maze = zeros(M,N,4);
      maze(:,:,1) = maze0(1:M,:);
      maze(:,:,2) = maze0(M+(1:M),:);
      maze(:,:,3) = maze0(2*M+(1:M),:);
      maze(:,:,4) = maze0(3*M+(1:M),:);
%     randomly open the maze at the bottom and top boundaries.      
      Xtop = randi([1,round(N/2)]);
      Xbottom = randi([round(N/2),N]);
      Xtop = 2;
      Xbottom = N-2;
      maze(1,[6 9],1)=0;
      maze(M,[2 7],3)=0;
      src = Xbottom -0.5i;
    end
    

%    mazeinfo = getmazeinfo(maze,M,N);
    mazeinfo = getconcisemazeinfo(maze,M,N);

    nsing  = mazeinfo.nsing;
    zsing  = mazeinfo.zsing;
    ncurve = mazeinfo.ncurve;
    ising  = mazeinfo.ising;
    clist  = mazeinfo.clist;
    iedge  = mazeinfo.iedge;
    iflag  = mazeinfo.iflag;

    cpars = cell(1,ncurve);
%   specify the boundary
    cparams = cell(1,ncurve);
    for icurve = 1:ncurve
      cparams{icurve}.ta = 0; % starting parameter value
      cparams{icurve}.tb = 1; % end parameter value
      cparams{icurve}.ifclosed = false; % whether each curve is closed or not
    end
    
%   starting and end points of each line segment
%    if isempty(ifplot)
%      figure; hold all; axis equal off
%    end

    for icurve = 1:ncurve
      cparams{icurve}.sing = ising(icurve,:);
      cpars{icurve}.v0=zsing(cparams{icurve}.sing(1));
      cpars{icurve}.v1=zsing(cparams{icurve}.sing(2));

%      if isempty(ifplot)
%	z1=zsing(cparams{icurve}.sing(1));
%	z2=zsing(cparams{icurve}.sing(2));
%	x=[real(z1) real(z2)];
%	y=[imag(z1) imag(z2)];
%	plot(x,y,'LineWidth',2.0);
%      end
    end

%    if isempty(ifplot)
%      xlim([0, N]);
%      ylim([0, M]);
%      hold off;
%      print -depsc maze.eps
%      ifplot=0;
%    end


    
%   structure contain information needed for singular points  
    nendpoint = 0;
    ncorner = 0;
    ntriple = 0;
    nquadruple = 0;

    singinfo = cell(1,nsing);
    for ising = 1:nsing
      singinfo{ising} = [];
      ifRcomp=1;
      isingcopy=0;

      nedge = iedge(ising);
      if nedge == 1
	nendpoint = nendpoint+1;
      elseif nedge == 2
	ncorner = ncorner+1;
      elseif nedge == 3
	ntriple = ntriple+1;
      elseif nedge == 4
	nquadruple = nquadruple+1;
      end
      
      curvelist = clist(ising,1:nedge);
      isstart = iflag(ising,1:nedge);
      
      singinfo{ising}.clist = curvelist;
      singinfo{ising}.isstart = isstart;
      singinfo{ising}.nedge = nedge;
      singinfo{ising}.ifRcomp = ifRcomp;
      singinfo{ising}.isingcopy = isingcopy;
    end

    [nendpoint ncorner ntriple nquadruple]
    
    osparams.cpars = cpars;
    osparams.cparams = cparams;
    osparams.ncurve = ncurve;
    osparams.nsing = nsing;
    osparams.singinfo = singinfo;
    osparams.zsing = zsing;
    osparams.src = src;
  end

  
  
  
  function osparams = ossetup1()
    osparams = [];
%   number of curves
    ncurve = 8; 
%   length of each line segment
    L=4; 
%   number of singular points
    nsing=9;
    
%   endpoints   
    z=zeros(nsing,1); 
    z(5)=0;
    z(6)=z(5)+L;
    z(3:4)=z(5)+exp(1i*2*pi/3*(1:-2:-1))*L;
    z(1:2)=z(3:4)-L;
    z(7:9)=z(6)+exp(1i*pi/3*(1:-1:-1))*L;

%   specify the boundary
    cparams = cell(1,ncurve);
    for icurve = 1:ncurve
      cparams{icurve}.ta = 0; % starting parameter value
      cparams{icurve}.tb = 1; % end parameter value
      cparams{icurve}.ifclosed = false; % whether each curve is closed or not
    end
    cparams{1}.sing=[1 3];
    cparams{2}.sing=[2 4];
    cparams{3}.sing=[3 5];
    cparams{4}.sing=[4 5];
    cparams{5}.sing=[5 6];
    cparams{6}.sing=[6 7];
    cparams{7}.sing=[6 8];
    cparams{8}.sing=[6 9];
    
    cpars = cell(1,ncurve);
%   starting and end points of each line segment  
    for icurve = 1:ncurve
      cpars{icurve}.v0=z(cparams{icurve}.sing(1));
      cpars{icurve}.v1=z(cparams{icurve}.sing(2));    
    end

%   structure contain information needed for singular points  
    singinfo = cell(1,nsing);
    for ising = 1:nsing
      singinfo{ising} = [];
      ifRcomp=1;
      isingcopy=0;
      if ising==1 % top left endpoint
%       curve ID starting or ending at the singular point
        clist = [1];
        isstart = [1];
%       number of edges      
        nedge = 1;
      elseif ising == 2 % bottom left endpoint
        clist = [2];
        isstart = [1];
        nedge = 1;
      elseif ising == 3 % top corner 
        clist = [1,3];
        isstart = [0,1];
        nedge = 2;
      elseif ising == 4 % bottom corner
        clist = [2,4];
        isstart = [0,1];
        nedge = 2;
      elseif ising == 5 % triple junction
        clist = [3,4,5];
        isstart = [0,0,1];
        nedge = 3;
      elseif ising == 6 % quadruple junction
        clist = [5,6,7,8];
        isstart = [0,1,1,1];
        nedge = 4;
      elseif ising == 7 % top right endpoint
        clist = [6];
        isstart = [0];
        nedge = 1;
      elseif ising == 8 % middle right endpoint
        clist = [7];
        isstart = [0];
        nedge = 1;
      elseif ising == 9 % bottom right endpoint
        clist = [8];
        isstart = [0];
        nedge = 1;
      end
  
      singinfo{ising}.clist = clist;
      singinfo{ising}.isstart = isstart;
      singinfo{ising}.nedge = nedge;
      singinfo{ising}.ifRcomp = ifRcomp;
      singinfo{ising}.isingcopy = isingcopy;
    end

    osparams.cpars = cpars;
    osparams.cparams = cparams;
    osparams.ncurve = ncurve;
    osparams.nsing = nsing;
    osparams.singinfo = singinfo;
    osparams.zsing = z;
  end

  
 
  function [z,zp,zpp] = funcurve(t,icurve,cpars,icase,curvetype)
%%funcurve
% return position, first and second derivatives of the curve
%
% Inputs:
% t - paramter values to evaluate these quantities
%
% Outputs:
% z - coordinates
% zp - first derivatives w.r.t. t
% zpp - second derivatives w.r.t. t
    if icase == 0 % line segment
      [z,zp,zpp] = linesegment(t,cpars);
    elseif icase == 1 % spiral
      [z,zp,zpp] = spiral(t,cpars);
    elseif icase == 2 % corners
      if curvetype == 1
	[z,zp,zpp] = linesegment(t,cpars);
      elseif curvetype == 2
	if mod(icurve,2) == 1
	  [z,zp,zpp] = parabola(t,cpars);
	else
	  [z,zp,zpp] = cosinecurve(t,cpars);
	end
      end
    elseif icase == 3 % branch points
      [z,zp,zpp] = linesegment(t,cpars);
    elseif icase == 4 % maze
      [z,zp,zpp] = linesegment(t,cpars);
    end
  end
  
  
  
  function [z,zp,zpp] = linesegment(t,cpars)
%%linesegment
% return position, first and second derivatives of the segment
% that passes through two endpoints v0,v1
%
% Inputs:
% t - paramter values [0,1] to evaluate these quantities
%
% Outputs:
% z - coordinates
% zp - first derivatives w.r.t. t
% zpp - second derivatives w.r.t. t

    z0 = cpars.v0;
    z1 = cpars.v1;

    islocal = -1;
    if isfield(cpars,'islocal')
      islocal = cpars.islocal; 
    end

    if islocal == -1
      z = z0 + (z1-z0)*t;
    else % local coordinates for RCIP
      z = (z1-z0)*t;
    end    
    zp = (z1-z0)*ones(size(t));
    zpp = zeros(size(t));
  end


  
  function [z,zp,zpp] = circulararc(t,cpars)
%%circulararc
% return position, first and second derivatives of a circular arc
% that passes through points v0 and v1 with angles theta0 and 
% theta1, respectively.
%
% Inputs:
% t - paramter values to evaluate these quantities
%
% Outputs:
% z - coordinates
% zp - first derivatives w.r.t. t
% zpp - second derivatives w.r.t. t

    r = cpars.r;
    zc = cpars.zc;
    theta0 = cpars.theta0;
    theta1 = cpars.theta1;

    islocal = -1;
    if isfield(cpars,'islocal')
      islocal = cpars.islocal;
    end

    dt = theta1-theta0;
    dt2 = dt^2;
  
    if islocal == -1
      theta = theta0+t*dt;

      z = r*exp(1i*theta);
      zp = 1i*dt*z;
      zpp = -dt2*z;

      z = zc+z;
    else
      if islocal == 0 % end point of the curve
        cx = r*cos(theta1);
        sx = r*sin(theta1);
      elseif islocal == 1 % starting point of the curve
        cx = r*cos(theta0);
        sx = r*sin(theta0);
      end
    
      st = sin(t*dt);
      stp = dt*cos(t*dt);
      stpp = -dt2*st;  
    
      ct = -2*sin(t*dt/2).^2;
      ctp = -dt*st;
      ctpp = -dt*stp;

      xs = -sx*st + cx*ct;
      ys =  cx*st + sx*ct;
      z = xs(:) + 1i*ys(:);
  
      xp = -sx*stp + cx*ctp;
      yp =  cx*stp + sx*ctp;
      zp = xp(:) + 1i*yp(:);
  
      xpp = -sx*stpp + cx*ctpp;
      ypp =  cx*stpp + sx*ctpp;
      zpp = xpp(:) + 1i*ypp(:);
    end
  end


  
  function [z,zp,zpp] = parabola(t,cpars)
%%parabola
% return position, first and second derivatives of a parabola
% defined by the formulas y=a(x-x0)^2+b, x=c1 t+c2;
%
% Inputs:
% t - paramter values to evaluate these quantities
%
% Outputs:
% z - coordinates
% zp - first derivatives w.r.t. t
% zpp - second derivatives w.r.t. t

    a = cpars.a;
    b = cpars.b;
    x0 = cpars.x0;
    c1 = cpars.c1;
    c2 = cpars.c2;

    islocal = -1;
    if isfield(cpars,'islocal')
      islocal = cpars.islocal;
    end
  
    if islocal == -1
      x = c1*t+c2;
      y = a*(x-x0).^2+b;
      
      xp = c1;
      yp = 2*a*(x-x0)*c1;
      
      xpp = zeros(length(t),1);
      ypp = 2*a*c1*c1;
      
    else
      if islocal == 0 % end point of the curve
        ctmp = c1+c2-x0;
      elseif islocal == 1 % starting point of the curve
        ctmp = c2-x0;
      end
    
      x = c1*t;
      y = a*c1*c1*t.^2+2*a*c1*ctmp*t;
      
      xp = c1;
      yp = 2*a*c1*c1*t+2*a*c1*ctmp;
      
      xpp = zeros(length(t),1);
      ypp = 2*a*c1*c1;
    end
    
    z = x+1i*y;
    zp = xp+1i*yp;
    zpp = xpp+1i*ypp;
  end

  
    
  function [z,zp,zpp] = spiral(tt,cpars)
%%spiral
% return position, first and second derivatives of a spiral
%
% Inputs:
% tt - paramter values to evaluate these quantities
%
% Outputs:
% z - coordinates
% zp - first derivatives 
% zpp - second derivatives 
%
%  z0 = cpars.z0; % starting point of the curve
    a = cpars.a; 
    b = cpars.b; % parameters a and b of the curve
    ta = cpars.ta; 
    tb = cpars.tb; % range of the parameter t

% islocal: -1 - regular points; 
%           1 - local coordinates close to the left endpoint;
%           0 - local coordinates close to the right endpoint.
    islocal = -1;
    if isfield(cpars,'islocal')
      islocal = cpars.islocal; 
    end
    if islocal == -1
      t=ta+tt;
      ebt=exp(1i*b*t);
%
      z=a*ebt;
      zp=1i*a*b*ebt;
      zpp=-a*b*b*ebt;
    elseif (islocal == 1) || (islocal == 0) 
% use stable and accurate forms to evaluate functions
% when points are close to one of the endpoints.
      t=tt;
      if islocal == 1
        t0=ta;
      elseif islocal == 0
        t0=tb;
      end
      e0=exp(1i*b*t0);
      st = sin(b*t);
      ctm1 = -2*sin(b*t/2).^2; % cos(t)-1
      t2 = tt+t0; % actual parameter value
      bt = b*t2;
      ebt=exp(1i*bt);
%
      z=a*e0*(ctm1+1i*st);
      zp=1i*a*b*ebt;
      zpp=-a*b*b*ebt;
    end    
  end



  
  function [z,zp,zpp] = cosinecurve(tt,cpars)
%% cosine curve
% return position, first and second derivatives of a cosine curve
%  z = z0+z1 (pi/2 t, cos(pi/2 t)), t \in [0,1]

    z0 = cpars.z0; z1 = cpars.z1;

    c = pi/2;
    
    islocal = -1;
    if isfield(cpars,'islocal')
      islocal = cpars.islocal; 
    end

    if islocal == -1
      xs = c*tt;
      ys = cos(c*tt);
      
      xp = c*ones(size(tt));
      yp = -c*sin(c*tt);

      xpp = zeros(size(tt));
      ypp = -c^2*ys;
    elseif islocal == 1 
      t = c*tt;
      ct = -2*sin(t/2).^2;

      xs = t;
      ys = ct;

      xp = c*ones(size(tt));
      yp = -c*sin(t);

      xpp = zeros(size(tt));
      ypp = -c^2*cos(t);
    elseif islocal == 0 
      t = c*tt;
      st = sin(t);
  
      xs = t;
      ys = -st;
      
      xp = c*ones(size(tt));
      yp = -c*cos(t);
      
      xpp = zeros(size(tt));
      ypp = -c^2*ys;
    end    
    

    z = xs(:) + 1i*ys(:);
    if islocal == -1
      z = z0+z1*z;
    else
      z = z1*z;
    end
    zp  = xp(:) + 1i*yp(:);   zp = z1*zp;
    zpp = xpp(:) + 1i*ypp(:); zpp = z1*zpp;
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Panelwise discretization subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [sinter,sinterdiff,npanfin]=panelinit12(nsub1,nsub2,npan)
%   divide [0,1] into npan+nsub1+nsub2 panels, with npan-2 coarse 
%   panels in the middle and nsub1/2+1 dyadically refined panels 
%   towards each endpoint
    npanfin=npan+nsub1+nsub2;
    sinter=zeros(npanfin+1,1);
    sinter(1:npan+1)=linspace(0,1,npan+1);
    sinterdiff=ones(npan+nsub1+nsub2,1)/npan;
%   *** left endpoint ***
    pnow=npan;
    for k=1:nsub1
      sinter(3:pnow+2)=sinter(2:pnow+1);
      sinter(2)=(sinter(1)+sinter(2))/2;   
      sinterdiff(3:pnow+1)=sinterdiff(2:pnow);
      sinterdiff(2)=sinterdiff(1)/2;
      sinterdiff(1)=sinterdiff(1)/2;   
      pnow=pnow+1;
    end
%   *** right endpoint ***
    for k=1:nsub2
      sinter(pnow+2)=sinter(pnow+1);
      sinter(pnow+1)=(sinter(pnow)+sinter(pnow+1))/2;   
      sinterdiff(pnow+1)=sinterdiff(pnow)/2;
      sinterdiff(pnow)=sinterdiff(pnow)/2;   
      pnow=pnow+1;    
    end
  end

  
  
  function [sinterfin,sinterdifffin,npanfin]=panelinit12b(nsub1,nsub2,...
    sinter,sinterdiff)
%   divide the first panel into nsub1 panels and the last panel into
%   nsub2 panels via dyadic refinement, keep the rest panels as is.
    npan=length(sinterdiff);
    npanfin=npan+nsub1+nsub2;
    sinterfin=zeros(npanfin+1,1);
    sinterfin(1:npan+1)=sinter;
    sinterdifffin=zeros(npan+nsub1+nsub2,1);
    sinterdifffin(1:npan)=sinterdiff;
%   *** left endpoint ***
    pnow=npan;
    for k=1:nsub1
      sinterfin(3:pnow+2)=sinterfin(2:pnow+1);
      sinterfin(2)=(sinterfin(1)+sinterfin(2))/2;   
      sinterdifffin(3:pnow+1)=sinterdifffin(2:pnow);
      sinterdifffin(2)=sinterdifffin(1)/2;
      sinterdifffin(1)=sinterdifffin(1)/2;   
      pnow=pnow+1;
    end
%   *** right endpoint ***
    for k=1:nsub2
      sinterfin(pnow+2)=sinterfin(pnow+1);
      sinterfin(pnow+1)=(sinterfin(pnow)+sinterfin(pnow+1))/2;   
      sinterdifffin(pnow+1)=sinterdifffin(pnow)/2;
      sinterdifffin(pnow)=sinterdifffin(pnow)/2;   
      pnow=pnow+1;    
    end
  end 

  
  
  function chnkr = chunkerfunc(fcurve,sinter,sinterdiff,T,W)
%CHUNKERFUNC create a chunker object for constructing the system matrix
% used in the forward recursion for computing the preconditioner R 
% in the RCIP method
% 
%
% Input: 
%   fcurve - function handle of the form
%               [z,zp,zpp] = fcurve(t)
%            where z, zp, zpp are size [dim,size(t)] arrays describing
%            position, first derivative, and second derivative of a curve
%            in dim dimensions parameterized by t.

    chnkr = []; % empty chunker

    npan = length(sinter)-1; % number of chunks
    ngl = length(T);

    np=ngl*npan;
    s=zeros(np,1);
    w=zeros(np,1);
    for k=1:npan
      myind=(k-1)*ngl+1:k*ngl;
      sdif=sinterdiff(k)/2;
      s(myind)=(sinter(k)+sinter(k+1))/2+sdif*T;
      w(myind)=W*sdif;
    end
    [z,zp,zpp]=fcurve(s) ;
    [zinter,~,~]=fcurve(sinter);
    nz=-1i*zp./abs(zp);
    wzp=w.*zp;
    awzp=w.*abs(zp);
    panlen=zeros(npan,1);
    for k=1:npan
      myind=(k-1)*ngl+1:k*ngl;
      panlen(k)=sum(awzp(myind));
    end
    
    chnkr.z=z;
    chnkr.zp=zp;
    chnkr.zpp=zpp;
    chnkr.s=s;
    chnkr.nz=nz;
    chnkr.w=w;
    chnkr.wzp=wzp;
    chnkr.awzp=awzp;
    chnkr.zinter=zinter;
    chnkr.sinter=sinter;
    chnkr.sinterdiff=sinterdiff;
    chnkr.npan=npan;
    chnkr.ngl=ngl;
    chnkr.panlen=panlen;
  end



  function chnkrtotal=mergechnkr(chnkr)
    ztot=[];nztot=[];wzptot=[];awzptot=[];zintertot=[];panlentot=[];
    stot=[];
    
    npan=0;
    for icurve=1:length(chnkr)
      ztot=[ztot; chnkr{icurve}.z];
      stot=[stot; chnkr{icurve}.s];
      nztot=[nztot; chnkr{icurve}.nz];
      wzptot=[wzptot;chnkr{icurve}.wzp];
      awzptot=[awzptot;chnkr{icurve}.awzp];
      npan0=chnkr{icurve}.npan;
      panlentot=[panlentot;chnkr{icurve}.panlen];
      npan=npan+npan0;
      zinter=chnkr{icurve}.zinter;
      zinter0=[zinter(1:npan0).';zinter(2:npan0+1).'];
      zintertot=[zintertot, zinter0];
    end
    chnkrtotal.z=ztot;
    chnkrtotal.s=stot;
    chnkrtotal.nz=nztot;
    chnkrtotal.wzp=wzptot;
    chnkrtotal.awzp=awzptot;
    chnkrtotal.npan=npan;
    chnkrtotal.ngl=chnkr{1}.ngl;
    chnkrtotal.panlen=panlentot;
    chnkrtotal.zinter=zintertot;
  end
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Construction of the sparse quadrature correction part of the system matrix
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [S,T] = osbuildmat_correction(chnkr,rpars)
% build the quadrature correction part of the system matrix for the open surface problem
% inputs: 
% rpars.k - array of wave numbers for each domain
% chnkr - array of chnkr objects specifying the whole boundary,
%         For constructing the preconditioner R in the RCIP method, chnkr
%         contains the type-b mesh on each edge
% outputs:
% S - system matrix correction for the single layer potential operator
% T - system matrix correction for the D' operator

    ncurve = length(chnkr);
    npanvec = zeros(ncurve,1);
    for i=1:ncurve
      npanvec(i)=chnkr{i}.npan;
    end
    
    omegaS = rpars.zkS;
    omegaT = rpars.zkT;
    ngl    = rpars.ngl;
    T16    = rpars.T;
    W16    = rpars.W;
    U      = rpars.U;
    ML0    = rpars.ML0;
    bctype = rpars.bctype;
    precond= rpars.precond;

%   total number of discretization points
    npan = sum(npanvec);
    np = sum(npanvec)*ngl;
    chnkrtotal = mergechnkr(chnkr);
    
%   now build the special quadrature correction part of the system matrix 
%    S = sparse(np,np);
%    T = sparse(np,np);
    Sdiag=zeros(ngl,ngl,npan);
    Ssupd=zeros(ngl,ngl,npan);
    Ssubd=zeros(ngl,ngl,npan);

    Tdiag=zeros(ngl,ngl,npan);
    Tsupd=zeros(ngl,ngl,npan);
    Tsubd=zeros(ngl,ngl,npan);

    for i=1:ncurve % curve id
%      indi = sum(npanvec(1:i-1))*ngl+(1:npanvec(i)*ngl);
      indi = sum(npanvec(1:i-1))+(1:npanvec(i));
%     block tridiagonal correction for each curve      
      [Sc, Tc] = osbuildmat_selfcorrection_sparse(chnkr{i},omegaS,omegaT, ...
					   T16,W16,U,ML0);
%      S(indi,indi)=Sc;
%      T(indi,indi)=Tc;
      Sdiag(:,:,indi)=Sc(:,:,:,1);
      Ssupd(:,:,indi)=Sc(:,:,:,2);
      Ssubd(:,:,indi)=Sc(:,:,:,3);

      Tdiag(:,:,indi)=Tc(:,:,:,1);
      Tsupd(:,:,indi)=Tc(:,:,:,2);
      Tsubd(:,:,indi)=Tc(:,:,:,3);
    end

%   now assemble all corrections in one sparse matrix for S and T.    
    Ssupd=Ssupd(:,:,2:end);
    Ssubd=Ssubd(:,:,1:end-1);
    
    Tsupd=Tsupd(:,:,2:end);
    Tsubd=Tsubd(:,:,1:end-1);

    u=Sdiag(:);
    u=[u;Ssubd(:)];
    u=[u;Ssupd(:)];

    v=Tdiag(:);
    v=[v;Tsubd(:)];
    v=[v;Tsupd(:)];
    npan=sum(npanvec);
%   generate the index arrays for sparse matrices.
%   first the main diagonal    
    [ind1,ind2,ind3]=ndgrid(0:ngl-1,0:ngl-1,0:npan-1);
    rind = 1+ind1(:)+ngl*ind3(:);
    cind = 1+ind2(:)+ngl*ind3(:);
%   sub-diagonal
    [ind1,ind2,ind3]=ndgrid(0:ngl-1,0:ngl-1,0:npan-2);
    rind = [rind;1+ngl+ind1(:)+ngl*ind3(:)];
    cind = [cind;1+ind2(:)+ngl*ind3(:)];
%   super-diagonal
    rind = [rind;1+ind1(:)+ngl*ind3(:)];
    cind = [cind;1+ngl+ind2(:)+ngl*ind3(:)];
%   build the sparse matrix in one call
    S = sparse(rind,cind,u,np,np);
    T = sparse(rind,cind,v,np,np);
  end




  function [Sc,Tc]=osbuildmat_selfcorrection_sparse(chnkr,omegaS,omegaT,T16,W16,U,ML0)
    z=chnkr.z;
    zp=chnkr.zp;
    nz=chnkr.nz;
    wzp=chnkr.wzp;
    awzp=chnkr.awzp;
    zinter=chnkr.zinter;
    sinterdiff=chnkr.sinterdiff;
    npan=chnkr.npan;
    ngl=length(T16);
    N=npan*ngl;
%   compute block tridiagonal correction matrices for log and hypersingular kernels  
    [Sc,Tc]=LogHypCinit_sparse(omegaS,omegaT,z,zp,nz,wzp,awzp,zinter, ...
			       sinterdiff,T16,W16,U,ML0,npan,N);
  end




  

  function [Sc,Tc]=LogHypCinit_sparse(omegaS,omegaT,z,zp,nz,wzp,awzp,zinter, ...
				      sinterdiff,T,W,A,LogC0,npan,N)
% *** Local panelwise evaluation ***
% *** iper=0,1 (0 is open arc, 1 is closed contour) ***
% *** isp=0,1 (0 is full, 1 is sparse) ***
% return correction matrices for log and hypersingular kernels.
% 
% Note well: 1. hypC is only for the operator T
%            2. uses Legendre polynomials with precomputed A instead of 
%               monomials.
    iper=0;
%    Sc=sparse(N,N);
%    Tc=sparse(N,N);
    ng=length(T);
    
    Sc=zeros(ng,ng,npan,3);
    Tc=zeros(ng,ng,npan,3);
    
    if iper==1
      kstart=1;
      kend=npan;  
    else
      kstart=2;
      kend=npan-1;  
    end
    
% compute U only once
    U=zeros(ng,ng,npan);
    for k=1:npan
      U(:,:,k)=wLCHS_src_precomp(zinter(k),zinter(k+1),z((k-1)*ng+(1:ng)));
    end
    
% *** central blocks ***
    for k=1:npan
      myind=(k-1)*ng+1:k*ng;
      dsd=sinterdiff(k)/2;
      S0=SoperC_targ(LogC0,omegaS,z(myind),z(myind),zp(myind),awzp(myind),dsd,1);
%      Sc(myind,myind)=S0;
      Sc(:,:,k,1)=S0;
      ztg=z(myind);
      nzt=nz(myind);
      [~,~,HypC]=wLCHS_target(zinter(k),zinter(k+1),ztg,z(myind), ...
			       nz(myind),wzp(myind),awzp(myind),U(:,:,k),1);
      T0=ToperC_targ(LogC0,HypC,omegaT,ztg,nzt,z(myind),zp(myind),nz(myind), ...
		     awzp(myind),dsd,1);
%      Tc(myind,myind)=T0;
      Tc(:,:,k,1)=T0;
    end
    
% *** superdiagonal blocks (targets to the left) ***
    for k=kstart:npan
      myinds=(k-1)*ng+1:k*ng;
      km1=mod(k-2,npan)+1;
      myindt=(km1-1)*ng+1:km1*ng;       

      alpha=sinterdiff(km1)/sinterdiff(k);
      [TMPL,accept,na]=WfrakLinit(A,-1-alpha,alpha,T);
      Ta=T+1;
      Tb=(1-T(accept))*alpha;
      LogC=TMPL./W'-log(abs(Ta'+Tb));

      dsd=sinterdiff(k)/2;
      S0=SoperC_targ(LogC,omegaS,z(myindt(accept)),z(myinds),zp(myinds), ...
		     awzp(myinds),dsd,0);
%      Sc(myindt(accept),myinds)=S0;
      Sc(accept,:,k,2)=S0;
      
      ztg=z(myindt);
      nzt=nz(myindt);   
      [~,~,HypC]=wLCHS_target(zinter(k),zinter(k+1),ztg,z(myinds), ...
			      nz(myinds),wzp(myinds),awzp(myinds),U(:,:,k),0);
      LogC1=zeros(ng,ng);
      LogC1(accept,:)=LogC;
      T0=ToperC_targ(LogC1,HypC,omegaT,ztg,nzt,z(myinds),zp(myinds), ...
		     nz(myinds),awzp(myinds),dsd,0);
      
%      Tc(myindt,myinds)=T0;
      Tc(:,:,k,2)=T0;
    end
    
% *** subdiagonal blocks (targets to the right) ***
    for k=1:kend
      myinds=(k-1)*ng+1:k*ng;
      kp1=mod(k,npan)+1;
      myindt=(kp1-1)*ng+1:kp1*ng;       
      alpha=sinterdiff(kp1)/sinterdiff(k);

      [TMPL,accept,na]=WfrakLinit(A,1+alpha,alpha,T);
      Ta=T-1;
      Tb=(T(accept)+1)*alpha;
      LogC=TMPL./W'-log(abs(Ta'-Tb));

      dsd=sinterdiff(k)/2;
      S0=SoperC_targ(LogC,omegaS,z(myindt(accept)),z(myinds),zp(myinds), ...
		     awzp(myinds),dsd,0);
%      Sc(myindt(accept),myinds)=S0;
      Sc(accept,:,k,3)=S0;
      
      ztg=z(myindt);
      nzt=nz(myindt);     
      [~,~,HypC]=wLCHS_target(zinter(k),zinter(k+1),ztg,z(myinds), ...
		  nz(myinds),wzp(myinds),awzp(myinds),U(:,:,k),0);
      LogC1=zeros(ng,ng);
      LogC1(accept,:)=LogC;
      T0=ToperC_targ(LogC1,HypC,omegaT,ztg,nzt,z(myinds),zp(myinds), ...
		     nz(myinds),awzp(myinds),dsd,0);
%      Tc(myindt,myinds)=T0;
      Tc(:,:,k,3)=T0;
    end
    
  end

  

  function S=SoperC_targ(LogC,omega,ztrg,z,zp,awzp,dsd,ifself)
% Build the correction matrix for the Helmholtz single layer potential.
    N=length(z);

    S = -besselj(0,omega*abs(ztrg-z.')).*LogC/pi;
    S = S.*awzp.';
    
    tmp=0.5i-(log(omega/2)+0.5772156649015328606)/pi;

    if ifself==1
      myind=1:N+1:N^2;
      S(myind)=(tmp-(LogC(myind)+log(dsd.*abs(zp)'))/pi).*awzp';
    end
  end

  
  
  function T=ToperC_targ(LogC,HypC,omega,ztrg,nzt,z,zp,nz,awzp,dsd,ifself)
% build the correction matrix for the Helmholtz T layer potential, i.e.,
% the normal derivative of the double layer potential.    

    Nt=length(ztrg);
    N=length(z);
    
    zdiff = ztrg-z.';
    TMP = omega*abs(zdiff);
    Ta1 = omega^2*besselj(1,TMP)./TMP.*real(nzt*nz');
    D1 = -real((ones(Nt,1)*nz.')./zdiff);
    D2 =  real((nzt*ones(1,N))./zdiff);
    Tb1 = TMP.^2.*besselj(2,TMP).*D1.*D2;
    T = -(Ta1+Tb1).*LogC/pi;
    T = T.*awzp.';

%   diagonal corrections
    if ifself==1
      gamma = 0.5772156649015328606;
      tmp = 0.5i-(log(omega/2)+gamma-0.5)/pi;

      myind=1:N+1:N^2;
      T(myind)=omega^2*0.5*(tmp-(LogC(myind)+log(dsd.*abs(zp)'))/pi).*awzp';
    end
    T = T-real(nzt.*HypC)/pi;
  end 


  

  function T=ToperC_near(LogC,HypC,omega,ztrg,nzt,z,nz,awzp)
% build the correction matrix for the Helmholtz T layer potential, i.e.,
% the normal derivative of the double layer potential.

%    nonself-interaction   

    Nt=length(ztrg);
    N=length(z);
    
    zdiff = ztrg-z.';
    TMP = omega*abs(zdiff);
    Ta1 = omega^2*besselj(1,TMP)./TMP.*real(nzt*nz');
    D1 = -real((ones(Nt,1)*nz.')./zdiff);
    D2 =  real((nzt*ones(1,N))./zdiff);
    Tb1 = TMP.^2.*besselj(2,TMP).*D1.*D2;

    T = -(Ta1+Tb1).*LogC/pi;
    T = T.*awzp.';
    T = T-real(nzt.*HypC)/pi;
  end 
  
  
  
  
  function [S,T]=osbuildmat_smooth_starind(chnkr,starind,omegaS,omegaT)
    z=chnkr.z(starind);
    nz=chnkr.nz(starind);
    awzp=chnkr.awzp(starind);
    N=length(z);
    
    zdiff=z-z.'; % target (column) - source (row)
    r=abs(zdiff);
    h0 = besselh(0,omegaT*r); 
    h1 = besselh(1,omegaT*r);

    zdiff = real(nz./zdiff);
    zdiff = zdiff.'.*zdiff;
    
    T = 0.5i*omegaT*(h1.*real(nz*nz')./r+r.*zdiff.*(2*h1-omegaT*r.*h0));
    T = T.*awzp.';
    T(1:N+1:end)=0;
    
    S = besselh(0,omegaS*r);
    S = 0.5i*S;
    S = S.*awzp.';
    S(1:N+1:end)=0;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Construction of the full system matrix for RCIP
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function M = osbuildmat(chnkr,rpars)
% build the system matrix for the open surface problem
% inputs: 
% rpars.k - array of wave numbers for each domain
% chnkr - array of chnkr objects specifying the whole boundary,
%         For constructing the preconditioner R in the RCIP method, chnkr
%         contains the type-b mesh on each edge
% outputs:
% S - system matrix for the single layer potential operator
% T - system matrix for the D' operator

    ncurve = length(chnkr);
    npanvec = zeros(ncurve,1);
    for i=1:ncurve
      npanvec(i)=chnkr{i}.npan;
    end
    
    omegaS = rpars.zkS;
    omegaT = rpars.zkT;
    ngl    = rpars.ngl;
    T16    = rpars.T;
    W16    = rpars.W;
    U      = rpars.U;
    ML0    = rpars.ML0;
    bctype = rpars.bctype;
    precond= rpars.precond;

    ifnearcorrection = false;
    if isfield(rpars,'ifnearcorrection')
      ifnearcorrection = rpars.ifnearcorrection;
    end
    
%   total number of discretization points
    np = sum(npanvec)*ngl;
    chnkrtotal = mergechnkr(chnkr);
    
%   now build the system matrix
    M = zeros(2*np);
%   smooth quadrature part
    ind1=1:np;
    ind2=(ind1)+np;
    [S,T]=osbuildmat_smooth(chnkrtotal,omegaS,omegaT);
    if bctype=='n'
      M(ind1,ind2) = T;
      M(ind2,ind1) = S; 
    else
      M(ind1,ind2) = S;
      M(ind2,ind1) = T;
    end
    
    for i=1:ncurve % curve id
      indi1 = sum(npanvec(1:i-1))*ngl+(1:npanvec(i)*ngl);
      indi2 = (indi1)+np;

      S0 = S(indi1,indi1);
      T0 = T(indi1,indi1);
%     diagonal block correction for each curve      
      [Sc, Tc] = osbuildmat_selfcorrection(S0,T0,chnkr{i},omegaS,omegaT, ...
					   T16,W16,U,ML0);

      M(indi1,indi1) =  M(indi1,indi1)-eye(size(Sc));
      if bctype=='n'
        M(indi1,indi2) =  M(indi1,indi2)+Tc;
        M(indi2,indi1) =  M(indi2,indi1)+Sc;
      else
        M(indi1,indi2) =  M(indi1,indi2)+Sc;
        M(indi2,indi1) =  M(indi2,indi1)+Tc;
      end
    end

%  near log/hypersingular corrections between end panels of two different
%  curves

    if ifnearcorrection
      for i=1:ncurve % target curve id
        indi1 = sum(npanvec(1:i-1))*ngl+(1:npanvec(i)*ngl);
        indi2 = (indi1)+np;
        for j=1:ncurve % source curve id
          if j~=i
            indj1 = sum(npanvec(1:j-1))*ngl+(1:npanvec(j)*ngl);
            indj2 = (indj1)+np;
            
            S0 = S(indi1,indj1);
            T0 = T(indi1,indj1);
	  
            [Sc,Tc] = osbuildmat_nearcorrection(S0,T0,chnkr{j},chnkr{i}, ...
						omegaS,omegaT,ngl);

            if bctype=='n'
              M(indi1,indj2) = M(indi1,indj2) + Tc;
              M(indi2,indj1) = M(indi2,indj1) + Sc;
            else
              M(indi1,indj2) = M(indi1,indj2) + Sc;
              M(indi2,indj1) = M(indi2,indj1) + Tc;
            end
          end
        end
      end
    end % end of if ifnearcorrection statement
    
  end



  
  function [S,T]=osbuildmat_smooth(chnkr,omegaS,omegaT)
    z=chnkr.z;
    nz=chnkr.nz;
    awzp=chnkr.awzp;
    N=length(z);
    
    zdiff = z-z.'; % target (column) - source (row)
    r = abs(zdiff);
    h0 = besselh(0,omegaT*r); 
    h1 = besselh(1,omegaT*r);
    zdiff = real(nz./zdiff);
    zdiff = zdiff.'.*zdiff;
    
%    T = 1i/2*omega*(h1.*nunup./r+omega*r2.*h2.*D1.*D2);
    T = 0.5i*omegaT*(h1.*real(nz*nz')./r+r.*zdiff.*(2*h1-omegaT*r.*h0));
    T = T.*awzp.';
    T(1:N+1:end) = 0;
    
    S = besselh(0,omegaS*r);
    S = 0.5i*S;
    S = S.*awzp.';
    S(1:N+1:end) = 0;
  end



  
  function [S,T]=osbuildmat_selfcorrection(S0,T0,chnkr,omegaS,omegaT,T16,W16,U,ML0)
    z=chnkr.z;
    zp=chnkr.zp;
    nz=chnkr.nz;
    wzp=chnkr.wzp;
    awzp=chnkr.awzp;
    zinter=chnkr.zinter;
    sinterdiff=chnkr.sinterdiff;
    npan=chnkr.npan;
    ngl=length(T16);
    N=npan*ngl;
%   compute block tridiagonal correction matrices for log and hypersingular kernels  
    [LogC,HypC]=LogHypCinit(z,nz,wzp,awzp,zinter,sinterdiff,T16, ...
			    W16,U,ML0,npan,N,0);
    S = SoperC(LogC,S0,chnkr,omegaS,ngl,N);
    T = ToperC(LogC,HypC,T0,chnkr,omegaT,ngl,N);
  end




  function [ML,M1]=LogHypCinit(z,nz,wzp,awzp,zinter,sinterdiff,T,W, ...
			       A,ML0,npan,N,isp)
% *** Local panelwise evaluation ***
% *** iper=0,1 (0 is open arc, 1 is closed contour) ***
% *** isp=0,1 (0 is full, 1 is sparse) ***
% return correction matrices for log and hypersingular kernels.
% 
% Note well: 1. hypC is only for the operator T
%            2. uses Legendre polynomials with precomputed A instead of 
%               monomials.
    iper=0;
    ML=zeros(N);
    M1=zeros(N);
    ng=length(T);

    if iper==1
      kstart=1;
      kend=npan;  
    else
      kstart=2;
      kend=npan-1;  
    end
    
% compute U only once
    U=zeros(ng,ng,npan);
    for k=1:npan
      U(:,:,k)=wLCHS_src_precomp(zinter(k),zinter(k+1),z((k-1)*ng+(1:ng)));
    end
    
% *** central blocks ***
    for k=1:npan
      myind=(k-1)*ng+1:k*ng;
      ML(myind,myind)=ML0;
      ztg=z(myind);
      nzt=nz(myind);
      [~,~,wcmpH]=wLCHS_target(zinter(k),zinter(k+1),ztg,z(myind), ...
	          nz(myind),wzp(myind),awzp(myind),U(:,:,k),1);
      M1(myind,myind)=nzt.*wcmpH;
    end
    
% *** superdiagonal blocks (targets to the left) ***
    for k=kstart:npan
      myinds=(k-1)*ng+1:k*ng;
      km1=mod(k-2,npan)+1;
      myindt=(km1-1)*ng+1:km1*ng;       

      alpha=sinterdiff(km1)/sinterdiff(k);
      [TMPL,accept,na]=WfrakLinit(A,-1-alpha,alpha,T);
      Ta=T+1;
      Tb=(1-T(accept))*alpha;
      MLloc=TMPL./W'-log(abs(Ta'+Tb));

      ML(myindt(accept),myinds)=MLloc;

      ztg=z(myindt);
      nzt=nz(myindt);   
      [~,~,wcmpH]=wLCHS_target(zinter(k),zinter(k+1),ztg,z(myinds), ...
		  nz(myinds),wzp(myinds),awzp(myinds),U(:,:,k),0);
      M1(myindt,myinds)=nzt.*wcmpH;
    end
    
% *** subdiagonal blocks (targets to the right) ***
    for k=1:kend
      myinds=(k-1)*ng+1:k*ng;
      kp1=mod(k,npan)+1;
      myindt=(kp1-1)*ng+1:kp1*ng;       
      alpha=sinterdiff(kp1)/sinterdiff(k);

      [TMPL,accept,na]=WfrakLinit(A,1+alpha,alpha,T);
      Ta=T-1;
      Tb=(T(accept)+1)*alpha;
      MLloc=TMPL./W'-log(abs(Ta'-Tb));

      ML(myindt(accept),myinds)=MLloc;

      ztg=z(myindt);
      nzt=nz(myindt);     
      [~,~,wcmpH]=wLCHS_target(zinter(k),zinter(k+1),ztg,z(myinds), ...
		  nz(myinds),wzp(myinds),awzp(myinds),U(:,:,k),0);
      M1(myindt,myinds)=nzt.*wcmpH;
    end
    
    M1=-M1/pi;
    if isp==1
      ML=sparse(ML);
      M1=sparse(M1);
    end
  end



  function S=SoperC(LogC,S0,chnkr,omega,ngl,N)
% Build the correction matrix for the Helmholtz single layer potential.
% Note a factor of 2 difference from the usual definition of 
% the Helmholtz fundamental solution. But "S" is standard thanks
% to an omitted factor of 2. At least according to Colton & Kress?

    z=chnkr.z;
    zp=chnkr.zp;
    nz=chnkr.nz;
    wzp=chnkr.wzp;
    awzp=chnkr.awzp;
    zinter=chnkr.zinter;
    sinterdiff=chnkr.sinterdiff;
    npan=chnkr.npan;
    if isreal(omega)
      S = -2*imag(S0).*LogC/pi;
    else
      S = -besselj(0,omega*abs(z-z.')).*LogC/pi;
      S = S.*awzp.';
    end
    
    tmp=0.5i-(log(omega/2)+0.5772156649015328606)/pi;

    dsd=sinterdiff(fix(((1:N)-1)/ngl)+1)'/2;
    myind=1:N+1:N^2;
    S(myind)=(tmp-(LogC(myind)+log(dsd.*abs(zp)'))/pi).*awzp';
  end



  function T=ToperC(LogC,HypC,T0,chnkr,omega,ngl,N)
% build the correction matrix for the Helmholtz T layer potential, i.e.,
% the normal derivative of the double layer potential.    
% Note a factor of 2 difference from the usual definition of
% the Helmholtz fundamental solution. But "T" is standard thanks
% to an omitted factor of 2. At least according to Colton & Kress?

    z=chnkr.z;
    zp=chnkr.zp;
    nz=chnkr.nz;
    wzp=chnkr.wzp;
    awzp=chnkr.awzp;
    zinter=chnkr.zinter;
    sinterdiff=chnkr.sinterdiff;
    npan=chnkr.npan;
    if isreal(omega)
      T=-2*imag(T0).*LogC/pi;
    else
      zdiff = z-z.';
      TMP = omega*abs(zdiff);
      Ta1 = omega^2*besselj(1,TMP)./TMP.*real(nz*nz');
      Dn = real((nz*ones(1,N))./zdiff);
      Tb1 = TMP.^2.*besselj(2,TMP).*Dn.'.*Dn;
      T = -(Ta1+Tb1).*LogC/pi;
      T = T.*awzp.';
    end
%   diagonal corrections  
    gamma = 0.5772156649015328606;
    tmp = 0.5i-(log(omega/2)+gamma-0.5)/pi;

    dsd=sinterdiff(fix(((1:N)-1)/ngl)+1)'/2;
    myind=1:N+1:N^2;
    T(myind)=omega^2*0.5*(tmp-(LogC(myind)+log(dsd.*abs(zp)'))/pi).*awzp';
    T = T+real(HypC);
  end 


  

  function [S,T]=osbuildmat_nearcorrection(S0,T0,chnkrsrc,chnkrtrg,omegaS,omegaT,ngl)
    npan    = chnkrsrc.npan;
    zsrc    = chnkrsrc.z;
    nzsrc   = chnkrsrc.nz;
    panlen  = chnkrsrc.panlen;
    zinter  = chnkrsrc.zinter;
    wzpsrc  = chnkrsrc.wzp;
    awzpsrc = chnkrsrc.awzp;
    
    zg   = chnkrtrg.z;
    nztrg  = chnkrtrg.nz;

    Ngtot  = length(zg);

    S = zeros(Ngtot,npan*ngl);
    T = zeros(Ngtot,npan*ngl);
    
    outlist=ones(Ngtot,1);

    panlen2=panlen.^2;
    dlim = 0.5;
%    dlim = 0.8;
    dlim2=dlim^2;
    
    panlist=(1:npan);
    targpan1=cell(1,npan);
    nspec=0;
    for k=1:Ngtot
      if outlist(k)
        zdiff=zg(k)-zsrc;
        xd=real(zdiff);
        yd=imag(zdiff);
        d2=xd.*xd+yd.*yd; % sqrt is expensive to compute, so avoid it.

        d2=reshape(d2,ngl,npan);
        d=min(d2)'./panlen2;
        tmp=d<dlim2;
        if nnz(tmp(panlist))>0
          panno=find(tmp)';
          for kk=panno
            targpan1{kk}=[targpan1{kk} k];
          end
          nspec=nspec+1;
        end
      end
    end
%    
    for kk=panlist
      k=targpan1{kk};
      if ~isempty(k)
        a=zinter(kk);
        b=zinter(kk+1);
        myind=(kk-1)*ngl+1:kk*ngl;
        nz=nzsrc(myind);
        zsc=zsrc(myind);
        awzp=awzpsrc(myind);
        wzp=wzpsrc(myind);

        U=wLCHS_src_precomp(a,b,zsc);
        
        ztg=zg(k);

        [LogC,CauC,HypC]=wLCHS_target(a,b,ztg,zsc,nz,wzp,awzp,U,0);
      	T(k,myind) = ToperC_near(LogC,HypC,omegaT,ztg,nztrg(k),zsc,nz,awzp);

        S(k,myind) = StargCloseCorrection(ztg,zsc,omegaS,LogC).*awzp';
      end
    end
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  RCIP subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function rcipstruct = rcipsetup(ngl,ndim,nedge,isstart,T,W) 
  % setup for the RCIP forward recursion
  % inputs:
  % ngl - number of Gauss-Legendre nodes
  % nedge - number of edges meeting at the corner
  % isstart - array of size nedge
  %           0 - means the edge ends at the corner w.r.t. the
  %           parameterization of the curve
  %           1 - means the edge starts from the corner w.r.t the
  %           parameterization of the curve
  % ndim - number of equations
  %
  % outputs:
  % Pbc - prolongation matrix, excluding trivial part
  % PbcL - prolongation matrix, including trivial part
  % PWbc - weighted prolongation matrix, excluding trivial part
  % starL, circL - bad and good indices for the local system matrix
  % starS, circS - bad and good indices for the preconditioner R
  %
    [IP,IPW]=IPinit(T,W);
    Pbc = Pbcinit2(IP,nedge,ndim);

    PbcL1=[];                         
    for i=1:nedge
      if isstart(i)
        PbcL1=blkdiag(PbcL1,IP,eye(ngl));
      else
        PbcL1=blkdiag(PbcL1,eye(ngl),IP);
      end
    end
    PbcL=PbcL1;
    for i=1:ndim-1
      PbcL=blkdiag(PbcL,PbcL1);
    end
    
    PWbc = Pbcinit2(IPW,nedge,ndim);
  % starL - bad indices for the system matrix M
  % circL - good indices for the system matrix M
    indg1 = 2*ngl + (1:ngl);
    indb1 = 1:2*ngl;
    indg0 = 1:ngl;
    indb0 = ngl + (1:2*ngl);
    starL0 = [];
    circL0 = [];  
    for i=1:nedge
      if isstart(i)
        starL0 = [starL0 indb1+3*(i-1)*ngl];
        circL0 = [circL0 indg1+3*(i-1)*ngl];
      else
        starL0 = [starL0 indb0+3*(i-1)*ngl];
        circL0 = [circL0 indg0+3*(i-1)*ngl];
      end
    end
    starL = starL0;
    circL = circL0;
    for i=1:ndim-1
      starL = [starL starL0+3*nedge*ngl]; % type-b mesh has 3 panels
      circL = [circL circL0+3*nedge*ngl];
    end

  % starS - bad indices for the preconditioner R
  % circS - good indices for the preconditioner R
    indg1 = ngl + (1:ngl);
    indb1 = 1:ngl;
    indg0 = 1:ngl;
    indb0 = ngl + (1:ngl);
    starS0 = [];
    circS0 = [];  
    for i=1:nedge
      if isstart(i)
        starS0 = [starS0 indb1+2*(i-1)*ngl];
        circS0 = [circS0 indg1+2*(i-1)*ngl];
      else
        starS0 = [starS0 indb0+2*(i-1)*ngl];
        circS0 = [circS0 indg0+2*(i-1)*ngl];
      end
    end
    starS = starS0;
    circS = circS0;
    for i=1:ndim-1
      starS = [starS starS0+2*nedge*ngl]; % type-c mesh has 2 panels
      circS = [circS circS0+2*nedge*ngl];
    end
    
    rcipstruct.Pbc=Pbc;
    rcipstruct.PbcL=PbcL;
    rcipstruct.PWbc=PWbc;
    rcipstruct.starL=starL;
    rcipstruct.starS=starS;
    rcipstruct.circL=circL;
    rcipstruct.circS=circS;
  end

  

  function [IP,IPW]=IPinit(T,W)
  % construct the prolongation matrix IP that maps function values
  % on n_{gl} Gauss-Legendre nodes on [-1,1] to function values at 
  % 2n_{gl} Gauss-Legendre, with shifted and scaled n_{gl}
  % Gauss-Legendre nodes on each subinterval [-1,0], [0,1], respectively.
  %
  % IPW is the weighted prolongation matrix acted on the left side. 
    ngl = length(T);
    A=ones(ngl);
    AA=ones(2*ngl,ngl);
    T2=[T-1;T+1]/2;
    W2=[W;W]/2;
    for k=2:ngl
      A(:,k)=A(:,k-1).*T;
      AA(:,k)=AA(:,k-1).*T2;   
    end
    IP=AA/A;
    IPW=IP.*(W2./W');
  end

  

  function Pbc=Pbcinit2(IP,nedge,ndim)
  % construct the nontrivial part of the prolongation matrix for the whole
  % system. Suppose that nedge is the number of edges meeting at the 
  % corner, ndim is the number of equations, n_{gl} is the number of 
  % Gauss-Legendre nodes on each chunk. Then Pbc is a block diagonal 
  % matrix with nedge diagonal blocks, where each block is of the size
  % 2n_{gl} \times n_{gl}

    Pbc=kron(eye(nedge*ndim),IP);
  end


  
  function [R,Rstor,Kstor]=Rcomp(nedge,ndim,nsub,rcipstruct,h01,isstart,...
				 fcurve,nse,rparslocal,T,W,ngl)
  % carry out the forward recursion for computing the preconditioner R 
  % in the RCIP method, using big block after big block convention.
  % This function also stores Rstor, and Kstor
  % for computing rhofin, i.e., the solution on the fine grid.
    starL = rcipstruct.starL;
    circL = rcipstruct.circL;
    starS = rcipstruct.starS;
    circS = rcipstruct.circS;
  
    Pbc  = rcipstruct.Pbc;
    PWbc = rcipstruct.PWbc;
  
  % size of the system matrix
    nsys = 3*ngl*nedge*ndim;
  
  % size of the preconditioner R
    nR = 2*ngl*nedge*ndim;

    Rstor = zeros(nR,nR,nse);
    Kstor = zeros(nsys,nsys,nse);
  
    sinter = zeros(4,nedge);
    sinterdiff = zeros(3,nedge);
  
    chnkrlocal=cell(1,nedge);
  
    for level=1:nsub
      h = h01(1,:)/2^(nsub-level);
      if level < nsub 
        for i=1:nedge
          if isstart(i)
            sinter(:,i) =  [0, 0.5, 1, 2]'*h(i);
            sinterdiff(:,i) =  [0.5, 0.5, 1]'*h(i);
          else
            sinter(:,i) = -[2, 1, 0.5, 0]'*h(i);
            sinterdiff(:,i) = [1, 0.5, 0.5]'*h(i);
          end
        end
      else 
  % at level=nsub, type c mesh contains two panels of lengths h0 and h1
  % type b mesh contains three panels of lengths h0/2, h0/2, and h1
        for i=1:nedge
          if isstart(i)
            sinter(:,i) =  [0, 0.5, 1, 1+h01(2,i)/h(i)]'*h(i);
            sinterdiff(:,i) =  [0.5*h(i); 0.5*h(i); h01(2,i)];
          else
            sinter(:,i) = -[1+h01(2,i)/h(i), 1, 0.5, 0]'*h(i);
            sinterdiff(:,i) = [h01(2,i); 0.5*h(i); 0.5*h(i)];
          end
        end 
      end
    % construct local chunks around the corner
      for i=1:nedge
        chnkrlocal{i} = chunkerfunc(fcurve{i},sinter(:,i),sinterdiff(:,i),T,W);
      end

    % construct the system matrix for local chunks
      if level==1
        rparslocal.ifnearcorrection = true;
      else
        rparslocal.ifnearcorrection = false;
      end
      MAT = osbuildmat(chnkrlocal,rparslocal);
      if level >= nsub-nse+1
        Kstor(:,:,nsub-level+1)=MAT;
      end
      MAT=MAT+eye(nsys);
      if level==1  
        R = myinv2a(MAT(starL,starL));
      end
      if level >= nsub-nse+1
        Rstor(:,:,nsub-level+1)=R;
      end
      R=SchurBana(Pbc,PWbc,MAT,R,starL,circL,starS,circS);
    end
  end



  function A=SchurBana(P,PW,K,A,starL,circL,starS,circS)
  % use matrix block inversion formula to recursively compute the
  % preconditioner R.
  % 
  % inputs: 
  % P, PW - nontrivial part (i.e., other than the identity matrices) 
  %         of prolongation and weighted prolongation matrices
  % K - the system matrix on a type-b mesh along each edge.
  %     the type-b mesh contains three chunks, with two small chunks
  %     close to the corner and one big chunk (twice of the small chunk
  %     in the parameter space) away from the corner.
  %
  % A - R_{i-1}, the (i-1)th iteration of R, R is on a type-c mesh.
  %     the type-c mesh contains two equal chunks, with the "bad" chunk
  %     close to the corner and the "good" chunk away from the corner.
  %
  %     R is basically a compression from the type-b mesh to the type-c
  %     mesh. When the recursion is done, R on the chunk closest to the
  %     corner contains all the information about the corner singularity in
  %     the sense that 1. during the solve stage, the system matrix is
  %     replace by R on that part of the curve; 2. during the eval stage,
  %     the transformed density rhotilde = R*rho can be integrated with any
  %     other smooth functions using the smooth quadrature rule, i.e., the
  %     Gauss-Legendre rule.
  %
  % starL, circL - indices of "bad" and "good" parts of the system matrix 
  %                bad part - two small chunks close to the corner
  %                good part - one big chunk away from the corner
  %
  % starS, circS - indices of bad and good parts of the preconditioner R
  %                bad part - the chunk close to the corner
  %                good part - the chunk away to the corner
  %
  % output:
  % A - R_{i}, the ith iteration of R.
  %
    VA=K(circL,starL)*A;
    PTA=PW'*A;
    PTAU=PTA*K(starL,circL);
%   stabilizing technique for this particular problem  
    DVAUI=myinv2(K(circL,circL)-VA*K(starL,circL));
%    *** improvement ***
%    DVAUI=myinv2b(K(circL,circL)-VA*K(starL,circL)); 
    DVAUIVAP=DVAUI*(VA*P);
    A(starS,starS)=PTA*P+PTAU*DVAUIVAP;
    A(circS,circS)=DVAUI;
    A(circS,starS)=-DVAUIVAP;
    A(starS,circS)=-PTAU*DVAUI;
  end


  
  function x=SchurBanaLs(K,A,rhs,starL,circL)
    AU=A*K(starL,circL);
    DVAU=K(circL,circL)-K(circL,starL)*AU;
    Arhs=A*rhs(starL);
    x=zeros(size(rhs));
    TMP=DVAU\(rhs(circL)-K(circL,starL)*Arhs);
    x(starL)=Arhs-AU*TMP;
    x(circL)=TMP;
  end

  
  
  function Mi=myinv2a(M)
    np=length(M)/2;
    B=M(1:np,np+1:2*np);
    C=M(np+1:2*np,1:np);
    Mi=[-inv(B*C) inv(C); inv(B) zeros(np)];
  end



  
  function Mi=myinv2b(M)
    np=length(M)/2;
    B=M(1:np,np+1:2*np);
    C=M(np+1:2*np,1:np);
    D=M(np+1:2*np,np+1:2*np);
    Mi=[ -C\D/B   inv(C); inv(B) zeros(np)];
  end

  
  
  function Mi=myinv2(M)
% compute the inverse of 2x2 block matrix
%
%  | A B |
%  | C D |
%
% when BC \approx -I, and A is close to singular

  np=length(M)/2;
  A=M(1:np,1:np);
  B=M(1:np,np+1:2*np);
  C=M(np+1:2*np,1:np);
  D=M(np+1:2*np,np+1:2*np);

  T=inv(A-B/D*C);

  Mi=[   T         -T*B/D;
      -D\C*T     inv(D)+D\C*T*B/D];
  end
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GMRES subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [x,it]=myGMRESR1(A,B,IG,R1,R2,R13,R4213,b,n,m,tol)
% *** GMRES with low-threshold stagnation control ***
    V=zeros(n,m+1);
    H=zeros(m);
    cs=zeros(m,1);
    sn=zeros(m,1);
    bnrm2=norm(b);
    V(:,1)=b/bnrm2;
    s=bnrm2*eye(m+1,1);
    for it = 1:m                                  
      it1=it+1;                                   
      w=matvec(A,B,IG,R1,R2,R13,R4213,V(:,it));
      for k = 1:it
        H(k,it)=V(:,k)'*w;
        w=w-H(k,it)*V(:,k);
      end
      H(it,it)=H(it,it)+1;
      wnrm2=norm(w);
      V(:,it1)=w/wnrm2;
      for k=1:it-1                                
        temp     = cs(k)*H(k,it)+sn(k)*H(k+1,it);
        H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
        H(k,it)  = temp;
      end
      [cs(it),sn(it)]=rotmatc(H(it,it),wnrm2);     
      H(it,it)= cs(it)*H(it,it)+sn(it)*wnrm2;
      s(it1) =-sn(it)*s(it);                      
      s(it)  = cs(it)*s(it);                         
      myerr=abs(s(it1))/bnrm2;
      if (myerr<=tol)||(it==m)                     
        disp(['predicted residual = ' num2str(myerr)])
        y=triu(H(1:it,1:it))\s(1:it);
        x=fliplr(V(:,1:it))*flipud(y);
        trueres=norm(x+matvec(A,B,IG,R1,R2,R13,R4213,x)-b)/bnrm2;
        disp(['true residual      = ',num2str(trueres)])
        break
      end
    end
  end


  
  function [x,it]=myGMRESR2(zkS,zkT,chnkr,S,T,Sc,Tc, ...
			    IG,R1,R2,R13,R4213,b,n,m,tol,bctype)
% *** GMRES with low-threshold stagnation control ***
    V=zeros(n,m+1);
    H=zeros(m);
    cs=zeros(m,1);
    sn=zeros(m,1);
    bnrm2=norm(b);
    V(:,1)=b/bnrm2;
    s=bnrm2*eye(m+1,1);
    for it = 1:m                                  
      it1=it+1;
      w=matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
		IG,R1,R2,R13,R4213,V(:,it),bctype);
      for k = 1:it
        H(k,it)=V(:,k)'*w;
        w=w-H(k,it)*V(:,k);
      end
      H(it,it)=H(it,it)+1;
      wnrm2=norm(w);
      V(:,it1)=w/wnrm2;
      for k=1:it-1                                
        temp     = cs(k)*H(k,it)+sn(k)*H(k+1,it);
        H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
        H(k,it)  = temp;
      end
      [cs(it),sn(it)]=rotmatc(H(it,it),wnrm2);     
      H(it,it)= cs(it)*H(it,it)+sn(it)*wnrm2;
      s(it1) =-sn(it)*s(it);                      
      s(it)  = cs(it)*s(it);                         
      myerr=abs(s(it1))/bnrm2;
      if mod(it,100) == 0
	[it myerr]
      end
      if (myerr<=tol)||(it==m)                     
        disp(['predicted residual = ' num2str(myerr)])
        y=triu(H(1:it,1:it))\s(1:it);
        x=fliplr(V(:,1:it))*flipud(y);
        trueres=norm(x+matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
			       IG,R1,R2,R13,R4213,x,bctype)-b)/bnrm2;
        disp(['true residual      = ',num2str(trueres)])
        break
      end
    end
  end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  GMRES with restart
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [x,iter] = myGMRESR3(zkS,zkT,chnkr,S,T,Sc,Tc, ...
                                     IG,R1,R2,R13,R4213,b,n, ...
				     restart,tol,bctype,maxit)
% *** GMRES with low-threshold stagnation control and restart option ***
    x = zeros(n, 1);               % Initial guess (usually zero)
    bnrm2 = norm(b);
    if bnrm2 == 0
      bnrm2 = 1.0;
    end
    r = b - x - matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
                        IG,R1,R2,R13,R4213,x,bctype);
    rnorm = norm(r);
    error = rnorm/bnrm2;
    if error < tol
      return
    end

    m0 = round(restart/2);
    m = restart+m0;    

    V = zeros(n,m+1);
    H = zeros(m);
    cs = zeros(m,1);
    sn = zeros(m,1);
    e1 = eye(m+1,1);
    s  = e1;
    
    pre = 0;
    for outerit = 1:maxit
      if outerit>1
	r = b - x - matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
                            IG,R1,R2,R13,R4213,x,bctype);
	rnorm = norm(r);
	pre = m0;
	for i=1:pre
	  V(:,1:pre)=V(:,restart+(1:pre));
	end
	H(1:pre,1:pre) = eye(pre);
	cs(1:pre) = 1;
	sn(1:pre) = 0;
	s = zeros(m+1,1);
      end
      
      V(:,1+pre) = r/rnorm;
      s(1+pre) = rnorm;
      
      for it = 1:restart                                  
        it1 = it + 1;
        w = matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
                    IG,R1,R2,R13,R4213,V(:,pre+it),bctype);
        for k = 1:it
          H(k+pre,it+pre) = V(:,k+pre)'*w;
          w = w - H(k+pre,it+pre)*V(:,k+pre);
        end
        H(pre+it,pre+it) = H(pre+it,pre+it) + 1;
        wnrm2 = norm(w);
        V(:,pre+it1) = w / wnrm2;
        for k = 1:it-1                                
          temp      =  cs(k+pre)*H(k+pre,it+pre) + sn(k+pre)*H(k+pre+1,it+pre);
          H(k+pre+1,it+pre) = -sn(k+pre)*H(k+pre,it+pre) + cs(k+pre)*H(k+pre+1,it+pre);
          H(k+pre,it+pre)   =  temp;
        end
        [cs(it+pre), sn(it+pre)] = rotmatc(H(it+pre,it+pre), wnrm2);     
        H(it+pre,it+pre) = cs(it+pre)*H(it+pre,it+pre) + sn(it+pre)*wnrm2;
        s(it1+pre) = -sn(it+pre)*s(it+pre);                      
        s(it+pre)  =  cs(it+pre)*s(it+pre);                         
        error = abs(s(it1+pre))/bnrm2;
        
	if error <= tol
          disp(['predicted residual = ', num2str(error)])
          y = triu(H(1:(it+pre),1:(it+pre)))\s(1:(it+pre));
          x = x + fliplr(V(:,1:(it+pre)))*flipud(y);
	  break;
	end
      end

      if error <= tol
	iter = it+(outerit-1)*m
	break;
      end
      y = triu(H(1:(restart+pre),1:(restart+pre)))\s(1:(restart+pre));
      x = x + fliplr(V(:,1:(restart+pre)))*flipud(y);
      r = b - x - matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
                          IG,R1,R2,R13,R4213,x,bctype);
      s(it1+pre) = norm(r);
      error = s(it1+pre)/bnrm2
      if error <= tol
	iter = it+(outerit-1)*m
	break;
      end
    end
  end




  function [x, flag, relres, iter] = gmres_restart_augmented(zkS,zkT,chnkr,S,T,Sc,Tc, ...
				     IG,R1,R2,R13,R4213,bctype, ...
				     b,tol,max_iter,restart,x0,keep)
    % GMRES with restart, keeping 'keep' vectors from the last iteration at each restart
    %
    % Inputs:
    %   A       - matrix or function handle representing the system matrix
    %   b       - right-hand side vector
    %   tol     - tolerance for convergence
    %   max_iter - maximum number of outer iterations (restarts)
    %   restart - number of inner iterations before restart
    %   x0      - initial guess (default: zero vector)
    %   keep    - number of vectors to retain from the previous restart cycle
    %
    % Outputs:
    %   x       - solution vector
    %   flag    - convergence flag (0 if converged, 1 otherwise)
    %   relres  - relative residual norm at the end
    %   iter    - number of iterations taken
    
    % Initialize parameters
    n = length(b);
    if nargin < 18, x0 = zeros(n, 1); end
    if nargin < 19, keep = 20; end  % Default to keeping the last 20 vectors if not specified
    x = x0;
    r = b - x - matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
                        IG,R1,R2,R13,R4213,x,bctype);
    beta = norm(r);
    relres = beta / norm(b);
    flag = 1;
    iter = 0;

    % Pre-allocate V and H with an extended size to retain "keep" vectors
    V = zeros(n, restart + keep + 1);
    H = zeros(restart + keep + 1, restart + keep);

    % Keep a record of last "keep" vectors across restarts
    last_V = [];

    % Restart loop
    while iter < max_iter && relres > tol
        % Initialize V and H for this cycle, including retained vectors
        if ~isempty(last_V)
            num_keep = min(keep, size(last_V, 2));
            V(:, 1:num_keep) = last_V(:, 1:num_keep);
        else
            num_keep = 0;
        end

        % Compute initial residual and start with orthonormal basis vector
	r = b - x - matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
                        IG,R1,R2,R13,R4213,x,bctype);
        beta = norm(r);
        V(:, num_keep + 1) = r / beta;
        
        % Inner loop (Arnoldi process)
        for j = 1:restart
            w = V(:,num_keep+j) + matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
                    IG,R1,R2,R13,R4213,V(:,num_keep+j),bctype);
            for i = 1:num_keep + j
                H(i, num_keep + j) = V(:, i)'*w;
                w = w - H(i, num_keep + j) * V(:, i);
            end
%            H(num_keep+j, num_keep+j) = H(num_keep+j, num_keep+j) + 1;
            H(num_keep + j + 1, num_keep + j) = norm(w);
            if H(num_keep + j + 1, num_keep + j) > 0
                V(:, num_keep + j + 1) = w / H(num_keep + j + 1, num_keep + j);
            end
            
            % Solve the least squares problem using QR factorization
            y = H(1:num_keep + j + 1, 1:num_keep + j) \ (beta * eye(num_keep + j + 1, 1));
            x = x0 + V(:, 1:num_keep + j) * y;
	    r = b - x - matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
				IG,R1,R2,R13,R4213,x,bctype);
	    relres = norm(r) / norm(b)
            iter = iter + 1;

            if relres < tol
                flag = 0;
                break;
            end
        end
        
        % Update x0 for the next restart
        x0 = x;

        % Store last "keep" vectors for the next restart
        last_V = V(:, max(1, num_keep + j - keep + 1):num_keep + j);

        if flag == 0
            break;
        end
    end

    if flag ~= 0 && relres <= tol
        flag = 0;
    end
  end
  



  
  function w=matvec(A,B,IG,R1,R2,R13,R4213,x)
    R1x=R1*x;
    BR1x=B*R1x;
    w=A*(R4213*BR1x+R2*x)-R13*BR1x-IG*R1x;
  end


  
  function w=matvec2(tol,zkS,zkT,chnkr,S,T,Sc,Tc, ...
		     IG,R1,R2,R13,R4213,x,bctype)
% S - the "bad" smooth part of the single layer potential system matrix
% it needs to be subtracted from the result of the FMM since the FMM
% computes the full smooth part of the system matrix

% Sc - the special quadrature correction part of the system matrix     
% it needs to added to the result of the FMM

% T - the "bad" smooth part of the D' system matrix    
% Tc - the special quadrature correction part of the system matrix
    R1x=R1*x;

%    BR1x=B*R1x; B=-S for the Neumann problem and -T for the Dirichlet problem
    if bctype == 'd'
      BR1x=T*R1x-Tc*R1x-Toper_fmm(tol,zkT,chnkr,R1x);
    elseif bctype == 'n'
      BR1x=S*R1x-Sc*R1x-Soper_fmm(tol,zkS,chnkr,R1x);
    end

    v=R4213*BR1x+R2*x;

    if bctype == 'd'
      Av=Soper_fmm(tol,zkS,chnkr,v)+Sc*v-S*v;
    elseif bctype == 'n'
      Av=Toper_fmm(tol,zkT,chnkr,v)+Tc*v-T*v;
    end

    w=Av-R13*BR1x-IG*R1x;
  end


  
  function [c,s]=rotmatc(a,b)
    if a==0
      c=0;
      s=1;
    else
      temp=b/a;
      c=1/sqrt(1+abs(temp)^2);
      s=temp*c;
    end
  end


  
  function [c,s]=rotmat(a,b)
    if b==0
      c=1;
      s=0;
    elseif abs(b)>abs(a)
      temp=a/b;
      s=1/sqrt(1+temp^2);
      c=temp*s;
    else
      temp=b/a;
      c=1/sqrt(1+temp^2);
      s=temp*c;
    end
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  kernel-split quadrature subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [WfrakL,accept,na]=WfrakLinit(A,trans,mscale,T)
% *** T is target vector, sources on canonical interval ***
% *** use Legendre polynomials instead of monomials!! ***
    ngl=length(T);
    accept=1:ngl;
    T=trans+mscale*T;
    accept=accept(abs(T)<2);
    na=length(accept);

    P=zeros(na,ngl+1);    
    Q=zeros(na,ngl);
    c0 = 2*(1:ngl-1)+1;
    c1 = 1./c0; % division is expensive, so do it only once.
    c2 = 1./(2:ngl);

    Ta=T(accept);
    upp=log(abs(1-Ta));
    loo=log(abs(1+Ta));
    P(:,1)=upp-loo;
    P(:,2)=2+Ta.*P(:,1);
    for k=1:ngl-1
      P(:,k+2)=(c0(k)*Ta.*P(:,k+1)-k*P(:,k))*c2(k); %c2(k)=1/(k+1),c0(k)=2*k+1
    end
    Q(:,1)=2*upp-(Ta+1).*P(:,1)-2;
    Q(:,2:ngl)=(P(:,1:ngl-1)-P(:,3:ngl+1)).*c1;
    WfrakL=Q*A; % matrix-vector product instead of inversion.
  end


  
  function U=wLCHS_src_precomp(a,b,zsc)
% computation about the source points, used before wLCHS_target
% O(p^3) operations
%
    Ng=length(zsc);
    cc=(b-a)/2;
    zsctr=(zsc-(b+a)/2)/cc; 
    U=legepols(zsctr,Ng-1);
    U=inv(U);
  end
  
  
  
  function [wcmpL,wcmpC,wcmpH,wcmpS]=wLCHS_target(a,b,ztg,zsc,nz,wzp,awzp,U,ifself)
% computation about the target points
% target ztg can now be a vector as well!!!
% O(p^2) operations per target for the Helmholtz kernels
%
% compute the correction to (nearly) log, singular, hypersingular, and supersingular
% integrals
% input arguments:
% a - left endpoint in complex variable form of the panel
% b - right endpoint in complex variable form of the panel
% ztg - target points in complex form
% zsc - source points in complex form
% nz - unit normal vectors at source points (in complex form)
% wzp - z'(t) w complex weights for the complexified contour integral
% awzp - |z'(t)|w real weights for the original line integral
% U - output of wLCHS_src_precomp
% ifself - whether the targets are on the source panel
%    
% output arguments: target - column, source - row
% wcmpL - correction for (nearly) logarithmically singualr integrals
% wcmpC - correction for (nearly) singular integrals
% wcmpH - correction for (nearly) hypersingular integrals
% wcmpS - correction for (nearly) supersingular integrals
%  
% *** zt is target vector in transformed plane ***
    Ng=length(zsc);
    nt=length(ztg);
    cc=(b-a)/2;
    zt=(ztg-(b+a)/2)/cc;zt=zt(:);

%   P_k is the Legendre polynomial of degree k    
    q=zeros(nt,Ng);   % q(k+1)=\int_{-1}^1 P_k(z)*log(z-zt) dz
    f=zeros(nt,Ng+1); % f(k+1)=\int_{-1}^1 P_k(z)/(z-zt) dz 

    if ifself==1
      upp=log(1-zt);
      loo=log(-1-zt);
      f(:,1)=upp-loo;
      ind=find(imag(f(:,1))<-1.1);
      f(ind,1)=f(ind,1)+1i*pi;
      ind=find(imag(f(1,:))>1.1);
      f(ind,1)=f(ind,1)-1i*pi;
    else
      ifleft = false(nt,1);
      zdiff = ztg-zsc.';
      xd = real(zdiff);
      yd = imag(zdiff);
      d2 = xd.*xd+yd.*yd;
      [~,ind]=min(d2,[],2);

      zsctr=(zsc-(b+a)/2)/cc;
      imzsctr=imag(zsctr);
      [~,myind]=max(abs(imzsctr));
      isconvex=sign(imzsctr(myind))==-1;
      for k=1:nt
        rezt=real(zt(k));
        imzt=imag(zt(k));
        costheta=xd(k,ind(k))*real(nz(ind(k)))+yd(k,ind(k))*imag(nz(ind(k)));
        if costheta<0
          ifleft(k)=true;                   % main rule
        end
        if abs(rezt)<1
          if isconvex                       % convex panel
            if costheta<0
              if imzt<0 && imzt>1.1*min(imzsctr)
                if imzt<mydepth(zsctr,rezt,Ng) 
                  ifleft(k)=false;          % convex panel exception
                end
              end
            end
          else                              % concave panel
            if costheta>0
              if imzt>0 && imzt<1.1*max(imzsctr)
                if imzt>mydepth(zsctr,rezt,Ng) 
                  ifleft(k)=true;           % concave panel exception
                end
              end
            end
          end
        end
      end
      
      gam=-1i*ones(nt,1);
      loggam=-0.5i*pi*ones(nt,1); % log(gam)
      gam(ifleft)=1i;
      loggam(ifleft)=0.5i*pi;
      f(:,1)=loggam+log((1-zt)./(gam.*(-1-zt)));
    end

    f(:,2)=2+zt.*f(:,1);
    for k=1:Ng-1
      f(:,k+2)=((2*k+1)*zt.*f(:,k+1)-k*f(:,k))/(k+1);
    end
    
%    q(:,1)=2*upp-(zt+1).*f(:,1)-2;
    P111=log(abs((1-zt).*(1+zt)));
    q(:,1)=P111-zt.*f(:,1)-2;
%   *** vectorization ***
    q(:,2:Ng)=(f(:,1:Ng-1)-f(:,3:Ng+1))./(3:2:2*Ng-1);
    
%   (nearly) log corrections
    tlog = log(abs((zsc.'-ztg)/cc));
    tlog(isinf(tlog))=0;
    wcmpL = imag((q*U)*cc.*conj(nz.'))./awzp.' - tlog;

%   (nearly) singular corrections    
    if nargout >= 2
      tcau = wzp.'./(zsc.'-ztg)/1i;
      tcau(isinf(tcau))=0;
      wcmpC = f(:,1:Ng)*U/1i - tcau;
    end

%   (nearly) hypersingular corrections    
    if nargout >= 3
      fp=zeros(nt,Ng);  % fp(k+1)=\int_{-1}^1 P'_k(z)/(z-zt) dz 
      fp(:,2)=f(:,1);
      for k=1:Ng-2
        fp(:,k+2) =(2*k+1)*f(:,k+1)+fp(:,k);
      end
%     Now add corrections so that fp(k+1)=\int_{-1}^1 P_k(z)/(z-zt)^2 dz
      c2=-(1./(1-zt)+(-1).^(0:Ng-1)./(1+zt));
      fp=fp+c2; 
      thyp = wzp.'./(zsc.'-ztg).^2/1i;
      thyp(isinf(thyp))=0;
      wcmpH = fp*U/(1i*cc) - thyp;
    end

%   (nearly) supersingular corrections        
    if nargout == 4
      fpp=zeros(nt,Ng); % fpp(k+1)=\int_{-1}^1 P"_k(z)/(z-zt) dz 
      for k=1:Ng-2
        fpp(k+2)=(2*k+1)*fp(k+1)+fpp(k);
      end
%     Now add corrections so that fpp(k+1)=\int_{-1}^1 P_k(z)/(z-zt)^3 dz
      c3=-0.5*(1./(1-zt).^2+(-1).^(1:Ng)/(1+zt).^2)...
        -0.25*(0:Ng-1).*(1:Ng).*(1./(1-zt)+(-1).^(1:Ng)./(1+zt));
      fpp=0.5*fpp+c3;
      tsup = wzp.'./(zsc.'-ztg).^3/1i;
      tsup(isinf(tsup))=0;
      wcmpS = fpp*U/1i/cc.^2 - tsup;
    end
  end



  function y=mydepth(zsctr,x,Ng)
    A=ones(Ng+2);
    xx=[-1;real(zsctr);1];
    for k=2:Ng+2
      A(:,k)=xx.*A(:,k-1);
    end
    coeff=A\[0;imag(zsctr);0];
    y=coeff(Ng+2);
    for k=Ng+1:-1:1
      y=y*x+coeff(k);
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  plot subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function fieldplot_abs(absa1,xg,yg,xylim,ngr,bctype,icase)
    F1=zeros(ngr);
    for k=1:ngr
      F1(1:ngr,k)=absa1((k-1)*ngr+(1:ngr));
    end
    set(0,'DefaultAxesFontSize',12)
    figure
    imagesc(xg,yg,F1);      
    colormap(jet)
    axis xy
    colorbar
    hold on

    xlabel('$x$','Interpreter','LaTeX','FontSize',17)      
    ylabel('$y$','Interpreter','LaTeX','FontSize',17)      

    axis equal
    axis(xylim)
    drawnow
    if icase == 1 && bctype == 'd'
      print -depsc spiralDfield.eps
    elseif icase == 1 && bctype == 'n'
      print -depsc spiralNfield.eps
    elseif icase == 2 && bctype == 'd'
      print -depsc cornerDfield.eps
    elseif icase == 2 && bctype == 'n'
      print -depsc cornerNfield.eps
    elseif icase == 3 && bctype == 'd'
      print -depsc branchDfield.eps
    elseif icase == 3 && bctype == 'n'
      print -depsc branchNfield.eps
    elseif icase == 4 && bctype == 'd'
      print -depsc mazeDfield.eps
    elseif icase == 4 && bctype == 'n'
      print -depsc mazeNfield.eps
    else
      print -depsc field.eps
    end
  end
  

  
  function fieldplot_real(u,xg,yg,xylim,ngr,bctype,icase)
    F1=zeros(ngr);
    for k=1:ngr
      F1(1:ngr,k)=u((k-1)*ngr+(1:ngr));
    end
    maxafield=max(max(abs(F1)))
    set(0,'DefaultAxesFontSize',12)
    figure(2)
    imagesc(xg,yg,F1,[-maxafield,maxafield]);      
    colormap(jet)
    axis xy
    colorbar
    hold on
    xlabel('$x$','Interpreter','LaTeX','FontSize',17)      
    ylabel('$y$','Interpreter','LaTeX','FontSize',17)      
%     if icase == 1
%       text(-8.6,-7.4,'(a)','FontSize',16)      
%     else
%       text(-9.6,-7.4,'(a)','FontSize',16)      
%     end
    axis equal
    axis(xylim)
    drawnow
    print -depsc field.eps
  end 
  
  
  
  function field_error_plot(udiff,xg,yg,xylim,ngr,bctype,icase)
    F1=zeros(ngr);
    for k=1:ngr
      F1(1:ngr,k)=abs(udiff((k-1)*ngr+(1:ngr)));
    end
    F1=log10(F1);
    maxa_err=max(max(F1));
    disp(' ')
%    disp(['original max field error = ',num2str(maxa_err)])
%    ind=find(F1>-9);
%    if ~isempty(ind)
%      nlarge=length(ind);
%      disp(['Warning: number of points with error > 1e-10 = ',num2str(nlarge)])
%      disp(['Perhaps something wrong in fieldcomp?'])
%    end
%    F1(ind) = log10(eps);
    
    maxa_err=max(max(F1));
    disp(['max field error = ',num2str(maxa_err)])
    set(0,'DefaultAxesFontSize',12)
    figure
    imagesc(xg,yg,F1,[log10(eps) maxa_err]);      
    colormap(flipud(pink(256)))
    axis xy
    colorbar
    hold on
    xlabel('$x$','Interpreter','LaTeX','FontSize',17)      
    ylabel('$y$','Interpreter','LaTeX','FontSize',17)
    axis equal
    axis(xylim)
    drawnow
    if icase == 1 && bctype == 'd'
      print -depsc spiralDfielderror.eps
    elseif icase == 1 && bctype == 'n'
      print -depsc spiralNfielderror.eps
    elseif icase == 2 && bctype == 'd'
      print -depsc cornerDfielderror.eps
    elseif icase == 2 && bctype == 'n'
      print -depsc cornerNfielderror.eps
    elseif icase == 3 && bctype == 'd'
      print -depsc branchDfielderror.eps
    elseif icase == 3 && bctype == 'n'
      print -depsc branchNfielderror.eps
    elseif icase == 4 && bctype == 'd'
      print -depsc mazeDfielderror.eps
    elseif icase == 4 && bctype == 'n'
      print -depsc mazeNfielderror.eps
    else
      print -depsc fielderror.eps
    end
  end


  
  function nsub_error_plot(myerr,nsubup,k,bctype,icase)  
    figure
    semilogy(1:k,myerr(1:k),'*')
    grid on
    axis([0 nsubup 1e-16 1])
    xlabel('number of subdivisions $n_{\rm sub}$','Interpreter','LaTeX', ...
           'FontSize',17)      
    ylabel('estimated relative error in  $u(r_{\rm targ})$','Interpreter', ...
           'LaTeX','FontSize',17)      
%    text(-16.6,2e-18,'(c)','FontSize',16)
    drawnow
    if icase == 0 && bctype == 'd'
      print -depsc lineDnsub.eps
    elseif icase == 0 && bctype == 'n'
      print -depsc lineNnsub.eps
    elseif icase == 1 && bctype == 'd'
      print -depsc spiralDnsub.eps
    elseif icase == 1 && bctype == 'n'
      print -depsc spiralNnsub.eps
    elseif icase == 2 && bctype == 'd'
      print -depsc cornerDnsub.eps
    elseif icase == 2 && bctype == 'n'
      print -depsc cornerNnsub.eps
    else
      print -depsc nsubstudy.eps
    end
  end


  function zk_iter_plot(iter,zk,bctype,icase)  
    n = length(zk);
    xmin = 40;
    if n == 7
      xmax = 4e3;
    elseif n == 8
      xmax = 7e3;
    elseif n == 9
      xmax = 1.4e4;
    elseif n == 10
      xmax = 2.7e4;
    end

    if icase == 0
      ymin = 1;
      ymax = 10;
    end

    if icase == 1
      ymin = 40;
      ymax = 100+10*(n-7);
    end

    figure
    loglog(zk,iter,'-o','LineWidth',2.0,'MarkerSize',10)
    grid on
    axis([xmin xmax ymin ymax])
    xlabel('ratio $\frac{L}{\lambda}$','Interpreter','LaTeX', ...
           'FontSize',17)      
    ylabel('number of iterations','Interpreter', ...
           'LaTeX','FontSize',17)      
%    text(-16.6,2e-18,'(c)','FontSize',16)
    drawnow
    if icase == 0 && bctype == 'd'
      print -depsc lineDiter.eps
    elseif icase == 0 && bctype == 'n'
      print -depsc lineNiter.eps
    elseif icase == 1 && bctype == 'd'
      print -depsc spiralDiter.eps
    elseif icase == 1 && bctype == 'n'
      print -depsc spiralNiter.eps
    elseif icase == 2 && bctype == 'd'
      print -depsc cornerDiter.eps
    elseif icase == 2 && bctype == 'n'
      print -depsc cornerNiter.eps
    elseif icase == 3 && bctype == 'd'
      print -depsc branchDiter.eps
    elseif icase == 3 && bctype == 'n'
      print -depsc branchNiter.eps
    end
  end


  function timing_plot(totaltime,gmrestime,buildtime,np,bctype,icase)  
    figure
    loglog(np,totaltime,'-o','LineWidth',2.0,'MarkerSize',10)
    hold on;
    loglog(np,gmrestime,'-d','LineWidth',2.0,'MarkerSize',10)
    loglog(np,buildtime,'-*','LineWidth',2.0,'MarkerSize',10)
%    loglog(np,evaltime,'-s','LineWidth',2.0,'MarkerSize',10)
    
    grid on
    legend('$T_{\rm total}$','$T_{\rm solve}$','$T_{\rm build}$',...
	   'Interpreter','LaTeX','FontSize',16,'Location','northwest')
%	   '$T_{\rm eval}$',...

    xlabel('$N$','Interpreter','LaTeX','FontSize',17)
    ylabel('Time (sec)','FontSize',17)
    drawnow
    if icase == 0 && bctype == 'd'
      print -depsc lineDtiming.eps
    elseif icase == 0 && bctype == 'n'
      print -depsc lineNtiming.eps
    elseif icase == 1 && bctype == 'd'
      print -depsc spiralDtiming.eps
    elseif icase == 1 && bctype == 'n'
      print -depsc spiralNtiming.eps
    elseif icase == 2 && bctype == 'd'
      print -depsc cornerDtiming.eps
    elseif icase == 2 && bctype == 'n'
      print -depsc cornerNtiming.eps
    elseif icase == 3 && bctype == 'd'
      print -depsc branchDtiming.eps
    elseif icase == 3 && bctype == 'n'
      print -depsc branchNtiming.eps
    end
  end


  function plot_boundary(ztot,xylim)
    set(0,'DefaultAxesFontSize',12)
    figure
    plot(real(ztot),imag(ztot),'b.')
    hold on
    axis equal
    axis(xylim)
    xlabel('$x$','Interpreter','LaTeX','FontSize',17)
    ylabel('$y$','Interpreter','LaTeX','FontSize',17)
    grid on
    drawnow
  end

  function plot_chunker(chnkr,xylim,icase)
    set(0,'DefaultAxesFontSize',12)
    c=colororder("gem12");
    lc=length(c);
    for i=1:length(chnkr)
      z=chnkr{i}.z;
%      plot(real(z),imag(z),'.','Color',c(mod(i-1,lc)+1,:),'LineWidth',2.0,'MarkerSize',10)
%      plot(real(z),imag(z),'.','Color',c(mod(i-1,lc)+1,:))
      plot(real(z),imag(z),'b.')
      hold on;
    end
    axis equal
    axis(xylim)
    xlabel('$x$','Interpreter','LaTeX','FontSize',17)
    ylabel('$y$','Interpreter','LaTeX','FontSize',17)
    grid on
    drawnow
    if icase == 0
      print -depsc line.eps
    elseif icase == 1
      print -depsc spiral.eps
    elseif icase == 2
      print -depsc corner.eps
    elseif icase == 3
      print -depsc branch.eps
    elseif icase == 4
      print -depsc maze.eps
    end
  end

  
  function plot_linesegments(cpars,xylim,icase)
    set(0,'DefaultAxesFontSize',12)
    c=colororder("gem12");
    lc=length(c);
    ic=0;
    for i=randperm(length(cpars))
      z0=cpars{i}.v0;
      z1=cpars{i}.v1;
      x=[real(z0), real(z1)];
      y=[imag(z0), imag(z1)];
%      plot(real(z),imag(z),'.','Color',c(mod(i-1,lc)+1,:),'LineWidth',2.0,'MarkerSize',10)
%      plot(x,y,'Color',c(mod(ic,lc)+1,:),'LineWidth',2.5)
      plot(x,y,'Color','blue','LineWidth',2.5)
      ic=ic+13;
      hold on;
    end
    axis equal
    axis(xylim)
    xlabel('$x$','Interpreter','LaTeX','FontSize',17)
    ylabel('$y$','Interpreter','LaTeX','FontSize',17)
    grid on
    drawnow
    if icase == 0
      print -depsc line.eps
    elseif icase == 1
      print -depsc spiral.eps
    elseif icase == 2
      print -depsc corner.eps
    elseif icase == 3
      print -depsc branch.eps
    elseif icase == 4
      print -depsc maze.eps
    end
  end

  
  function [arclen,zsrc,chnkr,osparams] = compute_curve_length(icase,helmosopts)
    osparams = ossetup(icase,helmosopts);
    curvetype = [];
    if isfield(osparams,'curvetype')
      curvetype = osparams.curvetype
    end
    ncurve = osparams.ncurve;
    cpars = osparams.cpars;
%   number of Gauss-Legendre nodes on each chunk
    ngl = 16;
    [T16,W16,~,~]=legeexps(ngl);
%   define functions for curves
    fcurve = cell(1,ncurve);
    for icurve=1:ncurve
      fcurve{icurve} = @(t) funcurve(t,icurve,cpars{icurve},icase,curvetype);
    end
    chnkr = cell(1,ncurve);
    npanvec = 100*ones(ncurve,1);
    arclen = zeros(ncurve,1);

    zsrc=[];
    for icurve = 1:ncurve
      [sinter,sinterdiff,~] = panelinit12(0,0,npanvec(icurve));
      chnkr{icurve} = chunkerfunc(fcurve{icurve},sinter,sinterdiff,T16,W16);
      arclen(icurve) = sum(chnkr{icurve}.awzp);
      zsrc=[zsrc;chnkr{icurve}.z];
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Legendre polynomial subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function p = legepols(x,n)
% return values of legendre polynomial p_0(x), p_1(x), ..., p_n(x)
    m=length(x);
    p=ones(m,n+1);
    if n == 0 
      return
    end
    p(:,2)=x;
    for k=1:n-1
      p(:,k+2)=( (2*k+1)*x.*p(:,k+1) - k*p(:,k) )/(k+1);
    end
  end


  
  function [nodes,whts,u,v] = legeexps(n)
% *** nodes and whts are Legendre nodes of weights of length n;
% *** u is n x n matrix converting function values at n Legendre
% *** nodes to Legendre expansion coefficients.
% *** v is the inverse of u
    [nodes, whts]=legewhts(n);

    if nargout < 3
      return
    end

    u=zeros(n);

    for i=1:n
      u(:,i)=legepols(nodes(i),n-1);
    end

    v=u';

    for i=1:n
      d=(2*i-1)/2;
      for j=1:n
        u(i,j)=v(j,i)*whts(j)*d;
      end
    end
  end

  
  
  function [ts, whts]=legewhts(n)
    eps=1d-15;
    h=pi/(2*n);

    ts=cos(h*(1:2:2*n-1));
    ts(floor(n/2)+1)=0;
    for i=1:floor(n/2)
      xk=ts(i);
      ifout=0;
      for k=1:10
        [pol,der,sum]=legepol_sum(xk,n);
        delta=-pol/der;
        xk=xk+delta;
        if abs(delta) < eps
          ifout=ifout+1;
        end
        if ifout == 3
          break;
        end
      end
      ts(i)=xk;
      ts(n-i+1)=-xk;
    end
    
    ts=sort(ts);
    ts=ts';

    if nargout == 1
      return
    end

    whts=zeros(n,1);
    for i=1:floor((n+1)/2)
      [pol,der,sum]=legepol_sum(ts(i),n);
      whts(i)=1/sum;
      whts(n-i+1)=whts(i);
    end
  end


  
  function [pol,der,sum]=legepol_sum(x,n)
    sum=0;
    pol=1;
    der=0;
    sum=sum+pol^2/2;
    if n==0
      return;
    end
  
    pol=x;
    der=1;
    sum=sum+pol^2*1.5;
    if n==1
      return
    end

    sum=0;

    pkm1=1;
    pk=x;
    sum=sum+pkm1^2/2;
    sum=sum+pk^2*1.5;

    pk=1;
    pkp1=x;
    for k=1:n-1
      pkm1=pk;
      pk=pkp1;
      pkp1=( (2*k+1)*x*pk-k*pkm1)/(k+1);
      sum=sum+pkp1^2*(k+1.5);
    end

    pol=pkp1;
    der=n*(x*pkp1-pk)/(x^2-1);
  end
 

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  maze related subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  function maze = generateMaze(M, N)
    % Initialize maze with walls
    % maze(y, x, direction): logical array indicating if wall exists (true) or not (false)
    maze = true(M, N, 4); % Directions: N=1, E=2, S=3, W=4

    % Initialize visited array
    visited = false(M, N);

    % Start carving from (1,1)
    [maze, ~] = carve_passages_from(1, 1, maze, visited);
%    drawMaze(maze);
  end

  function [maze, visited] = carve_passages_from(cx, cy, maze, visited)
    visited(cy, cx) = true; % Mark current cell as visited

    % Directions
    dx = [0, 1, 0, -1];
    dy = [-1, 0, 1, 0];
    opposite = [3, 4, 1, 2];

    directions = randperm(4); % Random order of directions

    for i = 1:4
      direction = directions(i);
      nx = cx + dx(direction);
      ny = cy + dy(direction);

      % Check if nx and ny are within bounds and not yet visited
      if ny >= 1 && ny <= size(maze, 1) && nx >= 1 && nx <= size(maze, 2) && ~visited(ny, nx)
        % Remove walls between current cell and next cell
        maze(cy, cx, direction) = false;
        maze(ny, nx, opposite(direction)) = false;

        % Recursive call
        [maze, visited] = carve_passages_from(nx, ny, maze, visited);
      end
    end
  end

  function drawMaze(maze)
    [M, N, ~] = size(maze);

    % Initialize figure
    figure; hold all; axis equal off;

    for y = 1:M
      for x = 1:N
        % Coordinates
        x1 = x - 1;
        x2 = x;
        y1 = M - y;
        y2 = M - y + 1;

        if maze(y, x, 1) % Top wall
          plot([x1, x2], [y2, y2],'LineWidth',2.0);
        end
        if maze(y, x, 2) % Right wall
          plot([x2, x2], [y1, y2],'LineWidth',2.0);
        end
        if maze(y, x, 3) % Bottom wall
          plot([x1, x2], [y1, y1],'LineWidth',2.0);
        end
        if maze(y, x, 4) % Left wall
          plot([x1, x1], [y1, y2],'LineWidth',2.0);
        end
      end
    end

    % Adjust axis limits
    xlim([0, N]);
    ylim([0, M]);
    hold off;
  end



  function mazeinfo = getmazeinfo(maze,M,N)
    nsing = (M+1)*(N+1)
    ncurve = sum(sum(maze(:,:,1)))+sum(sum(maze(:,:,2))) ...
	      +sum(maze(M,:,3))+sum(maze(:,1,4))
    
    ising = zeros(ncurve,2); % indices of endpoints for each curve
    clist=zeros(nsing,4); % indices of curves for each singular point
    iedge = zeros(nsing,1); % number of line segments connected to each singular point
%   flag showing whether the singular point is the starting point or the end point
    iflag = zeros(nsing,4); 
    
%   find the edges, direction always from left to right, or from south to north    
    icurve=0;
    for y=1:M
      for x=1:N
	% Coordinates
        x1 = x - 1;
        x2 = x;
        y1 = M - y;
        y2 = M - y + 1;
	
        if maze(y, x, 1) % Top wall
	  icurve=icurve+1;
	  inds=x1+1+y2*(N+1);
	  inde=x2+1+y2*(N+1);
	  ising(icurve,:) = [inds inde];

	  iedge(inds) = iedge(inds)+1;
	  clist(inds,iedge(inds))=icurve;
	  iflag(inds,iedge(inds))=1;

	  iedge(inde) = iedge(inde)+1;
	  clist(inde,iedge(inde))=icurve;
	  iflag(inde,iedge(inde))=0;
        end
	
        if maze(y, x, 2) % Right wall
	  icurve=icurve+1;
	  inds=x2+1+y1*(N+1);
	  inde=x2+1+y2*(N+1);
	  ising(icurve,:) = [inds inde];

	  iedge(inds) = iedge(inds)+1;
	  clist(inds,iedge(inds))=icurve;
	  iflag(inds,iedge(inds))=1;

	  iedge(inde) = iedge(inde)+1;
	  clist(inde,iedge(inde))=icurve;
	  iflag(inde,iedge(inde))=0;
        end

	if maze(y, x, 3) && y==M % Bottom wall
	  icurve=icurve+1;
	  inds=x1+1+y1*(N+1);
	  inde=x2+1+y1*(N+1);
	  ising(icurve,:) = [inds inde];

	  iedge(inds) = iedge(inds)+1;
	  clist(inds,iedge(inds))=icurve;
	  iflag(inds,iedge(inds))=1;

	  iedge(inde) = iedge(inde)+1;
	  clist(inde,iedge(inde))=icurve;
	  iflag(inde,iedge(inde))=0;
        end

	if maze(y, x, 4) && x==1 % Left wall
	  icurve=icurve+1;
	  inds=x1+1+y1*(N+1);
	  inde=x1+1+y2*(N+1);
	  ising(icurve,:) = [inds inde];    

	  iedge(inds) = iedge(inds)+1;
	  clist(inds,iedge(inds))=icurve;
	  iflag(inds,iedge(inds))=1;

	  iedge(inde) = iedge(inde)+1;
	  clist(inde,iedge(inde))=icurve;
	  iflag(inde,iedge(inde))=0;
        end
      end
    end
    
%   all grid points    
    [x,y]=ndgrid(0:N,0:M);
    z=x+1i*y;zsing=z(:);

    mazeinfo=[];
    mazeinfo.nsing  = nsing;
    mazeinfo.zsing  = zsing;
    mazeinfo.ncurve = ncurve;
    mazeinfo.ising  = ising;
    mazeinfo.clist  = clist;
    mazeinfo.iedge  = iedge;
    mazeinfo.iflag  = iflag;
  end

  
  
  
  function mazeinfo = getconcisemazeinfo(maze,M,N)
% remove artificial singular points and merge straight line segments    
    nsing = (M+1)*(N+1);
    ncurve = sum(sum(maze(:,:,1)))+sum(sum(maze(:,:,2))) ...
	     +sum(maze(M,:,3))+sum(maze(:,1,4));
    
    ising = zeros(ncurve,2); % indices of endpoints for each curve
    clist=zeros(nsing,4); % indices of curves for each singular point
    
%   find the edges, direction always from left to right, or from south to north    
    icurve=0;
    for y=1:M
      for x=1:N
	% Coordinates
        x1 = x - 1;
        x2 = x;
        y1 = M - y;
        y2 = M - y + 1;
	
	if maze(y, x, 4) && x==1 % Left wall
	  icurve=icurve+1;
	  inds=x1+1+y1*(N+1);
	  inde=x1+1+y2*(N+1);
	  ising(icurve,:) = [inds inde];    

	  clist(inds,3)=icurve;
	  clist(inde,4)=icurve;
        end

	if maze(y, x, 3) && y==M % Bottom wall
	  icurve=icurve+1;
	  inds=x1+1+y1*(N+1);
	  inde=x2+1+y1*(N+1);
	  ising(icurve,:) = [inds inde];

	  clist(inds,1)=icurve;
	  clist(inde,2)=icurve;
        end

        if maze(y, x, 1) % Top wall
	  icurve=icurve+1;
	  inds=x1+1+y2*(N+1);
	  inde=x2+1+y2*(N+1);
	  ising(icurve,:) = [inds inde];

	  clist(inds,1)=icurve;
	  clist(inde,2)=icurve;
        end
	
        if maze(y, x, 2) % Right wall
	  icurve=icurve+1;
	  inds=x2+1+y1*(N+1);
	  inde=x2+1+y2*(N+1);
	  ising(icurve,:) = [inds inde];

	  clist(inds,3)=icurve;
	  clist(inde,4)=icurve;
        end
      end
    end

    singremove = zeros(nsing,1);
    nh = sum(clist(:,1:2)>0,2);
    nv = sum(clist(:,3:4)>0,2);

    for i=1:nsing
      if nh(i) == 0 && nv(i) == 0
	singremove(i)=1;
      end
      if nh(i) == 2 && nv(i) == 0
	singremove(i)=1;
      end
      if nh(i) == 0 && nv(i) == 2
	singremove(i)=1;
      end
    end

    nsingnew = nsing-sum(singremove);

    curveremove = zeros(ncurve,1);
    for i=1:ncurve
      inds = ising(i,1);
      inde = ising(i,2);
      if singremove(inds) == 1 && singremove(inde) == 1
	curveremove(i) = 1; % throw away
      end
      if singremove(inds) == 1 && singremove(inde) == 0
	curveremove(i) = 2; % merge with 3, then throw away
      end
      if singremove(inds) == 0 && singremove(inde) == 1
	curveremove(i) = 3; % merge with 1 and 3, and change inde.
      end
    end

    ncurvenew = ncurve-sum(curveremove==1);

    for icurve=1:ncurve
      if curveremove(icurve) == 3
	inde = ising(icurve,2);
	currentcurve = icurve;
	while singremove(inde) == 1
	  if clist(inde,2) == currentcurve
	    nextcurve = clist(inde,1);
	  elseif clist(inde,4) == currentcurve
	    nextcurve = clist(inde,3);
	  end
	  assert(ising(nextcurve,1)==inde)
	  inde = ising(nextcurve,2);
	  currentcurve = nextcurve;
	end
	ising(icurve,2) = inde;
	if clist(inde,2) == currentcurve
	  clist(inde,2) = icurve;
	elseif clist(inde,4) == currentcurve
	  clist(inde,4) = icurve;
	end
	curveremove(icurve) = 0;
      end
    end

    ncurvenew = sum(curveremove==0);

    curveid = zeros(ncurve,1);
    jcurve=0;
    for icurve=1:ncurve
      if curveremove(icurve) == 0
	jcurve=jcurve+1;
	curveid(icurve)=jcurve;
      end
    end

    singid = zeros(nsing,1);
    jsing = 0;
    for i=1:nsing
      if singremove(i) == 0
	jsing=jsing+1;
	singid(i)=jsing;
      end
    end

%   all grid points    
    [x,y]=ndgrid(0:N,0:M);
    z=x+1i*y;z=z(:);

    clistnew = zeros(nsingnew,4);
    isingnew = zeros(ncurvenew,2);
    zsing = zeros(nsingnew,1);

    
    for i=1:ncurve
      if curveremove(i) == 0
	icurve=curveid(i);
	inds = ising(i,1);
	inde = ising(i,2);
	indsnew = singid(inds);
	indenew = singid(inde);

	isingnew(icurve,:) = [indsnew, indenew];
      end
    end

    for i=1:nsing
      if singremove(i) == 0
	inew = singid(i);
	for j=1:4
	   icurve = clist(i,j);
	   if icurve == 0
	     clistnew(inew,j) = 0;
	   else
	     icurvenew = curveid(icurve);
	     clistnew(inew,j) = icurvenew;
	   end
	end
	zsing(inew) = z(i);
      end
    end


    ising = isingnew;
    nsing = nsingnew;
    ncurve = ncurvenew;
    
    clist=zeros(nsing,4); % indices of curves for each singular point
    iedge = zeros(nsing,1); % number of line segments connected to each singular point
%   flag showing whether the singular point is the starting point or the end point
    iflag = zeros(nsing,4); 
    
    for icurve=1:ncurve
      inds = ising(icurve,1);
      inde = ising(icurve,2);

      iedge(inds) = iedge(inds)+1;
      clist(inds,iedge(inds))=icurve;
      iflag(inds,iedge(inds))=1;

      iedge(inde) = iedge(inde)+1;
      clist(inde,iedge(inde))=icurve;
      iflag(inde,iedge(inde))=0;
    end
      
%    figure; hold all; axis equal off;
%    for icurve=1:ncurvenew
%      z1 = zsing(ising(icurve,1));
%      z2 = zsing(ising(icurve,2));
%      plot([real(z1) real(z2)], [imag(z1) imag(z2)],'LineWidth',2.0);
%    end
%    hold off
    
    mazeinfo=[];
    mazeinfo.nsing  = nsing;
    mazeinfo.ncurve = ncurve;
    mazeinfo.ising  = ising;
    mazeinfo.zsing  = zsing;
    mazeinfo.clist  = clist;
    mazeinfo.iedge  = iedge;
    mazeinfo.iflag  = iflag;
  end

  
  
