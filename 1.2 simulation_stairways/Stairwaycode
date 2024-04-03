clear all;
close all;
tic



tol = 0.001;
maxits = 300;
dif = tol + 1000; % P_1 - P_0
numOfIterations = 0;


f_1 = 1/2; % distribution of cost
f_2 = 1/2;

theta_1 = 1; % marginal cost
theta_2 = 2;

d = 0.9; % discount rate per-stage


%% First Best Solution
qFB_1 = 1/(4*(theta_1)^2)
qFB_2 = 1/(4*(theta_2)^2)
vFB = 1/(1-d)*(f_1*(sqrt(qFB_1)-theta_1*qFB_1) + f_2*(sqrt(qFB_2)-theta_2*qFB_2))


if 1 % vgrid
    if 0 % uniform grid for v
        vmin = 0; % minimum, maximum vFB
        numv = 30; % grid points + 1
        % grid = (vFB-vmin)/vgrid; % grid
        vgrid = vmin:(vFB-vmin)/numv:vFB;
    end


    if 1   % nonuniform grid for v
        lower_vGrid = linspace(0, .4, 39);
        upper_vGrid = linspace(.401, vFB, 4);
        vgrid = [lower_vGrid upper_vGrid ];
    end

    vgrid = vgrid';
    Nv = size(vgrid,1);
end

%% Also Create the Grid for q, u

qmin=0;
numq=19;
qgrid= qmin:(qFB_1-qmin)/numq:qFB_1;
qgrid=qgrid';
[Nq,n]=size(qgrid);


umin=0;
umax=3;
numu=19;
ugrid= umin:(umax-umin)/numu:umax;
ugrid=ugrid';
[Nu,n]=size(ugrid);


%% Define a Function P(0), intialized as all zero function
p0 = zeros(Nv,1);

%% Update the P function throughout min
while dif > tol && numOfIterations < maxits
    for i = 1:Nv
        v0 = vgrid(i,1);
        p1cand = -999; % intialized a number for each data point; 
        % Do bruteforce checking; go through u_1,w_1,w_2,q_1,q_2
        for i_1 = 1:Nv
            w_1 = vgrid(i_1,1);
            pw_1=p0(i_1,1);
            for i_2 = 1:Nv
                w_2 = vgrid(i_2,1);
                pw_2 = p0(i_2,1);
                for j_1=1:Nq
                    q_1=qgrid(j_1,1);
                    for j_2=1:Nq
                        q_2=qgrid(j_2,1);
                        for l_1=1:Nu
                            u_1=ugrid(l_1,1);
                            u_2=2*v0-u_1-d*(w_1+w_2); %u_2 defined by (PK)
                            if u_1<0 ...
                                    || u_2<0 ...
                                    || u_1+d*w_1 < u_2+d*w_2+(theta_2-theta_1)*q_2 ...
                                    || u_2+d*w_2 < u_1+d*w_1+(theta_1-theta_2)*q_1
                                va = -99999999999;
                            else
                                va = f_1*(sqrt(q_1)-theta_1*q_1-u_1 + d*pw_1)...
                                    +f_2*(sqrt(q_2)-theta_2*q_2-u_2 + d*pw_2);
                            end
                            if va>p1cand
                                p1cand=va;
                                u1cand=u_1;
                                u2cand=u_2;
                                w1cand=w_1;
                                w2cand=w_2;
                                q1cand=q_1;
                                q2cand=q_2;
                            end
                        end
                    end
                end
            end
        end
        p1(i,1) = p1cand; %update the value function for each grid point i.e. optimal value for the programming per period)
        u1(i,1) = u1cand; %update the policy
        u2(i,1) = u2cand;
        w1(i,1) = w1cand;
        w2(i,1) = w2cand;
        q1(i,1) = q1cand;
        q2(i,1) = q2cand;
    end
    dif = norm(p1-p0)
    p0 = p1;
    numOfIterations = numOfIterations+1
end

if 1
    figure
    plot(vgrid, p1)
end

