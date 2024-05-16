%%% STL and robustness degree computing for multiple tasks (formulas)
%%% G[4,6](x1 >= 9) ^ F[10,12](x1 <= -10) ^ G[16,18](x1 <= -12 ^ F[22,24](x1 >= 13) ^ G[28,30](x1 <= -15)
%%%% Generation of random tasks. We pick one of those tasks and run SCP with initial point provided by the result in SQP15 
clc
clear
close all
rng("default")
load InitialX 
load InitialU 
x_init = xk;
u_init = uk';
m=1;
k_prime=2;
A=[0 1;-k_prime/m -1/5];
B=[0;1];
Q = [1 0;0 1]*0.01;
R = 0.001; %%% in R = 0.01 it's OK, however with 0.001 the signals are stabilized faster.
alpha = 1e-3; %%% Very sensitive to this parameter. greater than 1e-4, the task is not met. See and try with 0.001
xl = -20;
xu = 20;
ul = -40;
uu = 40;
% trust = 4;
x0 = [pi; -2];

N = 300;
x = sym('x',[2 N]);
k1 = 5;
dt = 0.1;
Nk = 5; %%% Number of iterations of the SCP
w_delta = 20; %%% Weight for penalization in the objective function due to trust region
Njtest = 10; %Number of tasks to test 
sigma_temp = 1.5;
sigma_x = 1.5;



%%% Temporal parameters for formulas %%% IT WAS CHANGED FOR THE EVENTUAL
%%% OPERATOR IN THE MIDDLE OF THE HORIZON: ALW EV ALW
l1_l = 4;
l1_u = 6;
l2_l = 10;
l2_u = 12;
l3_l = 16;
l3_u = 18;
l4_l = 22;
l4_u = 24;
l5_l = 28;
l5_u = 30;

%%% Spatial parameters for formulas
x1_limit1 = 9; 
x1_limit2 = -10; 
x1_limit3 = -12; 
x1_limit4 = 13; 
x1_limit5 = -15; 

  
x_total = [];
u_total = [];
rho_total = [];
Rho_degree_tasks = [];

tic

xk = x_init;
uk = u_init;
 
% xk(1,l1_l/dt:l1_u/dt) = x1_limit1*ones(1,l1_u/dt+1-l1_l/dt);
% xk(1,l2_l/dt:l2_u/dt) = x1_limit2*ones(1,l2_u/dt+1-l2_l/dt);

for j = 1:Njtest
    rho_sum1 = 0;
    rho_sum2 = 0;
    rho_sum3 = 0;
    rho_sum4 = 0;
    rho_sum5 = 0;
    
    flag = 0;
    while flag == 0
        l1_l_pert(j) = normrnd(l1_l/dt,sigma_temp/dt);
        l1_u_pert(j) = normrnd(l1_u/dt,sigma_temp/dt);
        if l1_l_pert(j) < l1_u_pert(j)
            flag = 1;
        end
    end

    flag = 0;
    while flag == 0
        l2_l_pert(j) = normrnd(l2_l/dt,sigma_temp/dt);
        l2_u_pert(j) = normrnd(l2_u/dt,sigma_temp/dt);
        if l2_l_pert(j) < l2_u_pert(j)
            flag = 1;
        end
    end

    flag = 0;
    while flag == 0
        l3_l_pert(j) = normrnd(l3_l/dt,sigma_temp/dt);
        l3_u_pert(j) = normrnd(l3_u/dt,sigma_temp/dt);
        if l3_l_pert(j) < l3_u_pert(j)
            flag = 1;
        end
    end

    flag = 0;
    while flag == 0
        l4_l_pert(j) = normrnd(l4_l/dt,sigma_temp/dt);
        l4_u_pert(j) = normrnd(l4_u/dt,sigma_temp/dt);
        if l4_l_pert(j) < l4_u_pert(j)
            flag = 1;
        end
    end

    flag = 0;
    while flag == 0
        l5_l_pert(j) = normrnd(l5_l/dt,sigma_temp/dt);
        l5_u_pert(j) = normrnd(l5_u/dt,sigma_temp/dt);
        if l5_l_pert(j) < l5_u_pert(j)
            flag = 1;
        end
    end

    x1_limit1_pert(j) = normrnd(x1_limit1,sigma_x);
    x1_limit2_pert(j) = normrnd(x1_limit2,sigma_x);
    x1_limit3_pert(j) = normrnd(x1_limit3,sigma_x);
    x1_limit4_pert(j) = normrnd(x1_limit4,sigma_x);
    x1_limit5_pert(j) = normrnd(x1_limit5,sigma_x);


    for i = 1:N

        if i >= l1_l_pert(j) && i <= l1_u_pert(j) %%% Specification: Signal X1 be always greater than 9 within l1_l to l1_u seconds
            rho_sum1 = rho_sum1 + exp(-k1*(x(1,i) - x1_limit1_pert(j)));
        end

        if i >= l2_l_pert(j) && i <= l2_u_pert(j) %%% Specification: Signal X1 eventually be less than -10 within l2_l and l2_u seconds
            rho_sum2 = rho_sum2 + exp(k1*(-x(1,i) + x1_limit2_pert(j)));
        end

        if i >= l3_l_pert(j) && i <= l3_u_pert(j) %%% Specification: Signal X1 be always greater than -12 within l3_l to l3_u seconds
            rho_sum3 = rho_sum3 + exp(-k1*(-x(1,i) + x1_limit3_pert(j)));
        end

        if i >= l4_l_pert(j) && i <= l4_u_pert(j) %%% Specification: Signal X1 eventually be greater than 13 within l4_l to l4_u seconds
            rho_sum4 = rho_sum4 + exp(k1*(x(1,i) - x1_limit4_pert(j)));
        end

        if i >= l5_l_pert(j) && i <= l5_u_pert(j) %%% Specification: Signal X1 be always less than -15 within l5_l to l5_u seconds
            rho_sum5 = rho_sum5 + exp(-k1*(-x(1,i) + x1_limit5_pert(j)));
        end

    end

    rho(j,1) = -1/k1*log(rho_sum1 + 1/rho_sum2 + rho_sum3 + 1/rho_sum4 + rho_sum5);
end


%%%%%% Clear these three lines below and set counter j in 1 again 
% clear
% load task7
% rho_total = [];
%%%%%%


for j = 1:Njtest
    for k =1:Nk
        %%% Inititalization

        rho_test = rho(j);


        %%%% Obtain the derivatives to compute the gradient

        for i = 1:N
            rho_grad(i) = diff(rho_test,x(1,i));
        end

        rho_orig = double(subs(rho_test,x(1,:),xk(1,:)));
        rho_grad_subs = double(subs(rho_grad,x(1,:),xk(1,:)));
        cvx_begin

        variable x(2,N)
        variable u(N)
        variable region
        expression rho_hat
        expression obj
        obj = 0;
        for i=1:N
            obj = obj + x(:,i)'*Q*x(:,i) + R*u(i)^2;
        end
        rho_hat = rho_orig + rho_grad_subs*(x(1,:)' - xk(1,:)');
        obj1 = rho_hat - alpha*obj - w_delta*region;
        maximize obj1
        subject to
        for i=1:N-1
            x(:,i+1) == x(:,i) + (A*x(:,i) + B*u(i))*dt;
        end
        x(:,1) == x0;
        ul <= u <= uu;
        xl <= x <= xu;
        region >= 0.8;
        for i=1:N
            norm(x(:,i) - xk(:,i)) <= region;
            norm(u(i) - uk(i)) <= region;
        end

        cvx_end

        xk = x;
        uk = u;

        x_total = [x_total;x];
        u_total = [u_total;u'];
        x = sym('x',[2 N]);
        Rho_degree = double(subs(rho_test,x(1,:),xk(1,:)));
        rho_total = [rho_total;Rho_degree];

        disp(['Task:',num2str(j),'  Iteration:',num2str(k),'  Robustness degree:',num2str(Rho_degree)])
        toc
    end
    

    Rho_degree_tasks = [Rho_degree_tasks,rho_total];
    
    namefile = strcat('Task: ',num2str(j));
    save(['Task',num2str(j),'.mat'])

    x_total = [];
    u_total = [];
    rho_total = [];
    xk = x_init;
    uk = u_init;

end

Rho_degree_mean = mean(Rho_degree_tasks,2);
figure;
plot(Rho_degree_tasks)
hold on
plot(Rho_degree_mean,'LineWidth',8)
xlabel('SCP iterations') 
ylabel('Robustness degree Rho') 
title(['Rho vs. Iterations for ' num2str(Njtest) ' tasks. ' 'STD: ' num2str(sigma_temp)])
grid on 

save AllTasks.mat

% figure; hold on
% t = linspace(0, N*dt,N);
% a1 = plot(t,x(1,:));
% M1 = "x1";
% a2 = plot(t, x(2,:));
% M2 = "x2";
% a3 = plot(t, u);
% M3 = "u";
% legend([a1,a2,a3], [M1, M2, M3]);
% grid on

time = toc;
toc











