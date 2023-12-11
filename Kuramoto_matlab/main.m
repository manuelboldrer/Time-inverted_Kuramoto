%/*-----------------------------------------------------
%|      Manuel Boldrer, PhD                            |
%|      Department of Cognitive Robotics               |
%|      Delft University of Technology                 |
%|                                                     |
%|      email: m.boldrer@tudelft.nl                    |
%|      https://manuelboldrer.github.io/               |                                                   
%-----------------------------------------------------*/
%The algorithms implemented in this code were used to generate
%part of the simulation results in the following papers:

%[1] Boldrer, M., Riz, F., Pasqualetti, F., Palopoli, L., & Fontanelli, D. (2021, December).
% Time-inverted Kuramoto dynamics for κ-clustered circle coverage.
% In 2021 60th IEEE Conference on Decision and Control (CDC) (pp. 1205-1211). IEEE.

%[2] Boldrer, M., Pasqualetti, F., Palopoli, L., & Fontanelli, D. (2022).
% Multiagent persistent monitoring via time-inverted kuramoto dynamics.
% IEEE Control Systems Letters, 6, 2798-2803.

%[3] Boldrer, M., Lyons, L., Palopoli, L., Fontanelli, D., & Ferranti, L. (2022).
% Time-Inverted Kuramoto Model Meets Lissajous Curves: 
% Multi-Robot Persistent Monitoring and Target Detection. 
% IEEE Robotics and Automation Letters, 8(1), 240-247.

clc, clear, close all

% Numerical simulation for the Kuramoto model
% theta_i'=Omega_i + K/N sum_j=1^N sin(theta_j-theta_i)

N=9;                         % Number of agents
p=5;                          % p determines the initial topology i.e., Ad
kk=1;                         % initialization
theta(:,1) = linspace(0,1,N); % agents' state
Ad = eye(N,N);                % Adjacency matrix

for m = 1:N  %Compute the Adjacency matrix depending on N and p.
    %notice that it can be any desired Ring topology
    TH = 1000;
    TH1 = 1000;
    TH2 = 0;
    TH3 = 1000;
    entrato = 0;
    entrato1 = 0;
    for n = 1:N
        if  m~=n && wrapTo360((theta(m,kk)-theta(n,kk))*180/pi) >= 0 &&...
                wrapTo360((theta(m,kk)-theta(n,kk))*180/pi) < TH % Rmaxt0(m)+Rmaxt0(n)
            TH = wrapTo360((theta(m,kk)-theta(n,kk))*180/pi) ;
            h = m;
            k = n;
            entrato1 = 1;
        end
    end

    if entrato1 == 1
        Ad(h,k) = 1;
        Ad(k,h) = 1;
    else
        for n = 1:N
            if m~=n && wrapTo360((theta(m,kk)-theta(n,kk))*180/pi) <= 0 &&...
                    abs(wrapTo360((theta(m,kk)-theta(n,kk))*180/pi)) < TH3 % Rmaxt0(m)+Rmaxt0(n)
                TH3 = abs(wrapTo360((theta(m,kk)-theta(n,kk))*180/pi));
                h = m;
                k = n;
            end
        end
        Ad(h,k) = 1;
        Ad(k,h) = 1;

    end
    Ad(h,k) = 1;
    Ad(k,h) = 1;

    for n = 1:N
        if m~=n && wrapTo360((theta(m,kk)-theta(n,kk))*180/pi) <= 0 &&...
                abs(wrapTo360((theta(m,kk)-theta(n,kk))*180/pi)) < TH1% Rmaxt0(m)+Rmaxt0(n)
            TH1 = abs(wrapTo360((theta(m,kk)-theta(n,kk))*180/pi));
            h1 = m;
            k1 = n;
            entrato = 1;
        end

    end
    if entrato == 1
        Ad(h1,k1) = 1;
        Ad(k1,h1) = 1;
    else
        for n = 1:N
            if m~=n && wrapTo360((theta(m,kk)-theta(n,kk))*180/pi) >= 0 &&...
                    abs(wrapTo360((theta(m,kk)-theta(n,kk))*180/pi)) > TH2% Rmaxt0(m)+Rmaxt0(n)
                TH2 = abs(wrapTo360((theta(m,kk)-theta(n,kk))*180/pi));
                h1 = m;
                k1 = n;
            end
        end
        Ad(h1,k1) = 1;
        Ad(k1,h1) = 1;

    end
end
%

Ad = Ad-eye(N); %Adjacency matrix

clearvars -except p N Ad
G  = graph(Ad);
Dd = zeros(N);
Dd = full(incidence(G));
L = full(laplacian(G));
neigh = cell(N,1);
iter = 2000;
diff = zeros(iter,N);
alpha = 4;
rs    = 0.45;
hend  = 10;
flag  = 1;
plotflag = 1;
figure

K = 0.1*N.*ones(N,1);
dt =  0.1;
for j=1:N
    theta(j,1) = (j-1)*2*pi*p/N+ 3*rand;
end
t = 0:dt:iter;

Rt = zeros(iter,1);
RR1 = zeros(N,iter);

%%
Omega = 0.05/alpha*ones(iter,N);%[0 0 0 0]';
for j=1:N
    thetaref(j,1) = (j-1)*2*pi*p/N + Omega(1,j).*dt;
end
ff = zeros(N,iter);
x= zeros(iter,N);
y= zeros(iter,N);

t = [0:dt:dt*iter]';

for kk = 1:iter
    for i = 1:N
        counte = 0;
        for j = 1:N
            if Ad(i,j) == 1
                counte = counte + 1;
                neigh{i}(counte,1) = j;
            end
        end

        diff(kk,i) = norm([cos(theta(i,kk))-cos(theta(neigh{i}(1),kk)),sin(theta(i,kk))-sin(theta(neigh{i}(1),kk))])+...
            norm([cos(theta(i,kk))-cos(theta(neigh{i}(2),kk)),sin(theta(i,kk))-sin(theta(neigh{i}(2),kk))]);
    end

    for j=1:N
        [x(kk,j),y(kk,j),z(kk,j),kappa(j)] =  fun(theta(j,kk),j,kk,[0,0],alpha);%2*cos(theta(j,kk));%1*cos(14*theta(j,kk))*cos(theta(j,kk));%+ 1*cos(theta(j,kk));%cos(theta(j,kk))^5;%1*(1-cos(theta(j,kk)))*cos(theta(j,kk));
        R2(j) = sqrt(x(kk,j)^2+y(kk,j)^2);%norm([ x(j),y(j)]);
    end
    countK = 0;
    for j = 1:N
        if Ad(1,j) == 1
            countK= countK+1 ;
            thetac(countK) = theta(j,kk);
        end
    end

    k1=kuramoto(theta(:,kk),K,N,Omega(kk,:)',Ad,Dd,thetaref);
    k2=kuramoto(theta(:,kk)+0.5*dt*k1,K,N,Omega(kk,:)',Ad,Dd,thetaref);
    k3=kuramoto(theta(:,kk)+0.5*dt*k2,K,N,Omega(kk,:)',Ad,Dd,thetaref);
    k4=kuramoto(theta(:,kk)+dt*k3,K,N,Omega(kk,:)',Ad,Dd,thetaref);
    theta(:,kk+1)=theta(:,kk)+(dt/6)*(k1+2*k2+2*k3+k4);

    ff(:,kk+1) = ((dt/6)*(k1+2*k2+2*k3+k4)); %4-th order Runge-Kutta method.

    R = 0;
    for j = 1:N
        R = R + (exp(1i*theta(j,kk)));
    end
    R1=zeros(N,1);
    for m = 1:N
        for j= 1:N
            if m==j || Ad(m,j) ==1
                R1(m) = R1(m) + exp(1i*theta(j,kk));
            end
        end
    end
    RR1(:,kk) = R1/N;
    RR(kk) = R/N;
    Rt(kk) = 1/N*norm(R);


    if kk ==1 %Plots
        co=1;
        for c=0:0.001:2*pi
            [xX(co),yY(co),zZ(co)] = fun(c,1,1,0,1);
            co=co+1;
        end
        figure(1)
        axis equal
        axis off
        plot3((xX),(yY),(zZ),'r','linewidth',1)

    end
    hold on

%%
    if plotflag == 1 %Plots

        if 1
            ag = scatter3(x(kk,:),y(kk,:),z(kk,:),'filled','MarkerFaceColor',DodgerBlue);
            for j = 1:N

                neigh1 = [0,0];
                cnt = 0;
                for jj = 1:N
                    if Ad(j,jj) == 1
                        cnt   = cnt+1;
                        neigh1(cnt) = theta(jj,kk);
                    end
                end
                co = 1;

            end
            hold on
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'ztick',[])
           
            [xo,yo] = gplot(Ad,[x(kk,:)',y(kk,:)']);
            edges = plot(xo,yo,'color',DarkOrange,'linewidth',1);




            axis equal
            set(gcf,'color','w')
        end
        box on
        drawnow
                F(kk) = getframe(gcf); %#ok<*SAGROW>

        delete(edges)
        if kk <iter
            delete(ag)
        end
    end
end


%%
% % %

video = VideoWriter('TIK6.avi','Motion JPEG AVI');
    video.Quality = 100;
    video.FrameRate = 1/(dt/4);
    open(video)
    writeVideo(video,F(1:1525))
    close(video)





