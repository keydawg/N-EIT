% $$$ function recon3(Rat_number,EIT_serial)

files=dir('Data/*.mat');
files={files.name};

for ff=1:length(files)

load('Jacobian_Tet_GSK_P9.mat', ...
     'BV0','J');

load('prt_ring_all.mat');

  load(['Data/' files{ff}]);
 
    interval=[5000:5400];
    
 
    [~,~,sel] = intersect(Prt_0,prt_0_all,'rows','stable');
    
    Noise = 1e-6*dZ(4800:5030,:);
    dZ  = 1e-6*dZ(interval,:);
    J=J(sel,:);
    BV0=BV0(sel);
    
 plot(dZ)
 drawnow;
    dZ(:,BV0<0) = -dZ(:,BV0<0);
    Noise(:,BV0<0) = -Noise(:,BV0<0);
    
   % plot(dZ)

    %%
%      figure;
%      scatter(abs(BV0),max(abs(dZ)));
%      xl=xlim; yl=ylim;
%      hold on
%      plot([0,max([xl,yl])],[0,max([xl,yl])],'-k')
%      drawnow;
%     % $$$ title(sprintf('Rat %03i EIT %s',Rat_number,EIT_serial))
    % $$$ clear BV BV0
    %%

    n_T = size(dZ,1);
    n_N = size(Noise,1);
    n_J = size(J,1);

    x=dZ;
    y=Noise;

    %%

    tic

    [U,S,V] = svd(J,'econ');
    disp(sprintf('SVD done: %.2f min.',toc/60))

    lambda = logspace(-16,-2,4000);

    [X,cv_error] = tikhonov_CV(J,x',lambda,n_J,U,S,V);
    disp(sprintf('X done: %.2f min.',toc/60))

    UtNoise = U'*y';
    sv = diag(S);

    % $$$ [~,opt] = min(cv_error);
    opt = zeros(size(cv_error,2),1);
    for i=1:size(cv_error,2)
        e = cv_error(:,i);
        f = (e(1:end-2)>=e(2:end-1))&(e(3:end)>e(2:end-1));
        if any(f)
            opt(i) = find(f,1,'last')+1;
        else
            [m,opt(i)] = min(e);
        end
    end

    SD = zeros(size(X));
    for i=1:size(SD,2)
        sv_i = sv+lambda(opt(i))./sv;
        SD(:,i) = std(V*(diag(1./sv_i)*UtNoise),0,2);
        disp(num2str([i,toc/60]))
    end


    save(['recon_no_bigs' files{ff}(1:end-5) '.mat'],...
         'X','cv_error','lambda','SD','interval','-v7.3');
end
