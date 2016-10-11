
function MD_coll_int_3
    global lambda_v lambda N cae Ree Ret gamma lambdaw
    gamma = 0.5; %surface coverage
    cae = 100; % <1 then attraction is dominate
    Ree = 1; % <1 then viscous forces dominate
    Ret = 100; % <1 then viscous forces dominate
    Len = 30; %length of square
    lambda = 1/2;
    lambdaw=1;
    
    %initial R 
    [R_wall, N_wall] = initialization_wall(Len,Len); % N_wall= total wall particle
    N = round(num_particle(Len,Len,N_wall) - N_wall);
    gamma = N*(pi * lambda^2)/(Len*Len);
    %R = initialization(N,Len,Len,lambda);
    R = initialization_lattice(Len,Len);
    R = [R; R_wall];

    %eqilibrium setup
    lambda_v = 0; % shear = 0 when lambda_v = 0;
    disp('equilibriation');
    equil_steps=100;
    [R_equilb,equil_end]=equilibrium(Len,R,equil_steps);

    % stepper for position update
    disp('simulation');
    lambda_v = 1; %  shear = 0 when lambda_v = 0;
    steps = 5; 
    delta_t=10e-5;
    R_new = stepper(Len,R_equilb{equil_end},steps,delta_t);
    
    %visualization
    title1=['Equilib Init Cae =',num2str(cae),' and Ree=',num2str(Ree)];
    temp=R_equilb{1};
    scatter_vis(1,temp(:,1),temp(:,2),Len,Len,title1);

    title2=['Equilib End Cae =',num2str(cae),' and Ree=',num2str(Ree)];
    temp=R_equilb{equil_end};
    scatter_vis(2,temp(:,1),temp(:,2),Len,Len,title2);
  
    temp=R_new{2};
    title3='Step = 1';
    scatter_vis(3,temp(:,1),temp(:,2),Len,Len,title3);
   
    temp=R_new{steps};
    title4=(['Step =', num2str(steps)]);
    scatter_vis(4,temp(:,1),temp(:,2),Len,Len,title4);
     
    
    movie_writer(R_new,steps,5,Len,Len);
    save('R_Equilb.mat', 'R_equilb');
    save('R_new.mat', 'R_new');
   
end  

function [p,n_wall] = initialization_wall(x_bound,y_bound)
% [x,y,phi,r (lambda), theta];
global lambda

p(:,1) = linspace(-x_bound/2,x_bound/2,x_bound/lambda);
p(:,2) = y_bound/2;
ln = length(p);

for i = 1:ln
    r0=3;
    p(i+ln,1) = p(i,1);
    p(i+ln,2) = -p(i,2);
    p(i,3) = 4*pi*rand - 2*pi; %Phi
    p(i+ln,3) = p(i,3);
    p(i,4) = 3/r0;%r (lambda)
    p(i+ln,4) = p(i,4);
    p(i,5) = 0;%theta
    p(i+ln,5) = p(i,5);
end
n_wall = length(p);
end

function [n] = num_particle(x_bound,y_bound,nwall)
global gamma lambda
n1 = gamma * x_bound * y_bound / (pi * lambda^2);
n=ceil(sqrt(n1))*ceil(sqrt(n1))+nwall;
end

function [p] = initialization(n,x_bound,y_bound,lambda)
% creates a population of particles in x and y bound with qualties of
% [x,y,phi,r (lambda), theta];
p = zeros(n,5);
for i = 1:n
    r0=3;
    p(i,1) = rand*(-x_bound) + x_bound/2;%x
    p(i,2) = rand*(-y_bound) + y_bound/2;%
    p(i,3) = 4*pi*rand - 2*pi; %Phi
    p(i,4) = lambda;%r
    p(i,5) = 0;%theta
end
end

function [p] = initialization_lattice(x_bound,y_bound)
    % creates a population of particles in x and y bound with qualties of
    % [x,y,phi,r (lambda), theta];

    global lambda N
    p = zeros(N,5);
    xpos=linspace(-x_bound/2+lambda,x_bound/2-lambda,ceil(sqrt(N)));
    ypos=linspace(-y_bound/2+2*lambda,y_bound/2-2*lambda,ceil(sqrt(N)));
    [x,y] = meshgrid(xpos,ypos);
    x(1:2:end,:)=x(1:2:end,:)+(x_bound/(2*ceil(sqrt(N))));
    coords = [reshape(x,1,numel(x)); reshape(y,1,numel(y))];
    for i=1:1:N
        p(i,1) = coords(1,i);
        p(i,2) = coords(2,i);
    end

    for i = 1:N
        p(i,3) = 4*pi*rand - 2*pi; %Phi
        p(i,4) = lambda;%r
        p(i,5) = 0;%theta
    end
end

function [F, T] = forces3(p)
global N cae

F = zeros(length(p),2);
T = zeros(length(p),1);

for i = 1:N
    for j = i+1:length(p)
        if j ~= i
                                  
            L = sqrt((p(j,1)-p(i,1))^2 + (p(j,2)-p(i,2))^2);
                        
            if L <= (p(i,4)+p(j,4))
                A = (-1/(L^11)+1);
                Fx = A * (p(j,1)-p(i,1)) / L; %Hard particle potential
                Fy = A * (p(j,2)-p(i,2)) / L;
                Txy = 0;
            else
                B = (p(i,4)*p(j,4))^2 * ((1/cae) * cos(2*(p(i,3)+p(j,3)))...
                    /(L^5) - ((1 + cos(p(i,5))) * (1 + cos(p(j,5)))) / (L^4));
                
                Fx = B * (p(j,1)-p(i,1)) / L;
           
                Fy = B * (p(j,2)-p(i,2)) / L;
                
                Txy = (p(i,4)*p(j,4))^2 * sin(2*(p(i,3)+p(j,3)))/ (L^4);
            end
                           
            F(i,1) = Fx + F(i,1);
            F(i,2) = Fy + F(i,2);
            F(j,1) = -Fx + F(j,1);
            F(j,2) = -Fy + F(j,2);
            T(i,1) = Txy + T(i,1);
            T(j,1) = -Txy + T(j,1);
                      
            end                         
        end       
end
F = F(1:N,:);% WHY ISNT THIS IN THE IMAGE FORCE CODE?
T = T(1:N,:);  
end

function [Fi, Ti] = forces_image(p,Len)
global N cae

p = p(1:N,:);
Fi = zeros(N,2);
Ti = zeros(N,1);
pi = zeros(N,5);

for i = 1:N
    for j = 1:N
        if j ~= i                 
                      
            pi(j,:) = image(p(j,:),Len);
            Li = sqrt((pi(j,1)-p(i,1))^2 + (pi(j,2)-p(i,2))^2);
            
            if pi(j,1) == 0   % ignoring irrelevent images
                Fix = 0;
                Fiy = 0;
                Tixy = 0;
            else
                if Li <= (p(i,4)+pi(j,4))
                    A = (-1/(Li^11)+1);
                    Fix = A * (pi(j,1)-p(i,1)) / Li; %Hard particle potential
                    Fiy = A * (pi(j,2)-p(i,2)) / Li;
                    Tixy = 0;
                else
                    B = (p(i,4)*pi(j,4))^2 * ((1/cae) * cos(2*(p(i,3)+pi(j,3)))...
                        /(Li^5) - (1 + cos(p(i,5))) * (1 + cos(pi(j,5))) / (Li^4));

                    Fix = B * (pi(j,1)-p(i,1)) / Li;

                    Fiy = B * (pi(j,2)-p(i,2)) / Li;

                    Tixy = (p(i,4)*pi(j,4))^2 * sin(2*(p(i,3)+pi(j,3)))/ (Li^4);
                end

            end
                           
            Fi(i,1) = Fix + Fi(i,1);
            Fi(i,2) = Fiy + Fi(i,2);
            Ti(i,1) = Tixy + Ti(i,1);
                        
        end       
    end
end
end

function [Results] = stepper(Len,R,steps,delta_t)
     Results=cell(steps+1,1);
     Results{1}=R; 
    for i = 1:1:steps

        [F, T] = forces3(R);
        [Fi, Ti] = forces_image(R,Len);
        F_final = F+Fi;
        T_final = T+Ti;
        v = velocity(F_final,R);   
        w=angular_velocity(T_final);
       

        for j = 1:1:length(R)
            R(j,1) = R(j,1) + v(j,1)*delta_t;
            R(j,2) = R(j,2) + v(j,2)*delta_t;
            R(j,1) = periodic_BC(R(j,1),Len);  
        end

        Results{i+1}=R;% recording all positions at every time step after
        %equilibriation.

        vis=rem(i,10);
        if vis==0
           disp(int2str(i));
        end
    end
end

function [Re,i] = equilibrium(Len,R,steps)
    %initial equilibriation    
    Re=cell(steps+1,1);
    Re{1}=R; 
    
    [F, T] = forces3(R);
    [Fi, Ti] = forces_image(R,Len);
    F_final = F+Fi;
    T_final = T+Ti;
    v = velocity(F_final,R);  
    w=angular_velocity(T_final);
    v_max=mean(max(v,[],1));
    delta_t=0.1/v_max;
    i=1;
    
    while delta_t < 5*10e-5
        [F, T] = forces3(R);
        [Fi, Ti] = forces_image(R,Len);
        F_final = F+Fi;
        T_final = T+Ti;
        v = velocity(F_final,R);   
        w=angular_velocity(T_final);
        v_max=max(max(v,[],1));
        delta_t=0.1/v_max;
        for j = 1:length(R)
            R(j,1) = R(j,1) + v(j,1)*delta_t;
            R(j,2) = R(j,2) + v(j,2)*delta_t;
            R(j,1) = periodic_BC(R(j,1),Len);  
        end
        
        Re{i+1}=R;% recording all positions at every time during
        %equilibriation.

        vis=rem(i,10);
        if vis==0
           disp(int2str(i));
           disp(delta_t);
        end
       i=i+1;    
       if i==steps
           delta_t=1;
       end
    end
end

function [Ri_element] = image(R_element,Len)

Ri_element = R_element;

if R_element(1,1) > 0; % image within [-L, L].. box within [-L/2, L/2]
    Ri_element(1,1) = R_element(1,1)-Len;
elseif R_element(1,1) < 0
    Ri_element(1,1) = R_element(1,1)+Len;
else
    Ri_element(1,1) = 0;
end

end

function [Rx] = periodic_BC(Rx,Len)

Rx = Rx - Len*round(Rx/Len); % periodic BC

end

function [v] = velocity(F,R)
global lambda_v Ree N lambdaw
v = zeros(length(R),2);

for i = 1:N
    v(i,1) = Ree*F(i,1) + lambda_v*lambdaw*R(i,2);
    v(i,2) = Ree*F(i,2);
end
for i = N+1:length(R) % for wall particle velocity
    v(i,1) = lambda_v*lambdaw*R(i,2);
    v(i,2) = 0;
end
end

function [w] = angular_velocity(T)
global N Ret
w = zeros(N,1);

for i = 1:N
    w(i,1) = Ret*T(i,1) -1;
end
end

function [E] = Energy(p)
global N cae

E = zeros(N,1);

for i = 1:N
    for j = i+1:N
        if j ~= i
                                  
            L = sqrt((p(j,1)-p(i,1))^2 + (p(j,2)-p(i,2))^2);
                        
            if L <= (p(i,4)+p(j,4))
                Ei = (1/(L^10)); %Hard particle potential
                                
            else
                Ei = (p(i,4)*p(j,4))^2 * ((-1/4)*(1/cae) * cos(2*(p(i,3)+p(j,3)))...
                    /(L^4) + (1/3)*(1 + cos(p(i,5))) * (1 + cos(p(j,5))) / (L^3));
                                
            end
            
            E(i,1) = Ei + E(i,1);
            E(j,1) = -Ei + E(j,1);
            
        end                         
    end       
end
end

function [] = scatter_vis(num,x,y,lenx,leny,name)
    figure(num)
    hold on
    scatter(x,y,10);
    grid on
    xlim([-2*lenx/2 2*lenx/2])
    ylim([-2*leny/2 2*leny/2])
    title(name);
    xlabel('L/r0')
    ylabel('L/r0')  
end

function [] = movie_writer(pos_cell,num_steps, fig_num, lenx, leny)
   myObj=VideoWriter('results','MPEG-4')
   open(myObj);
    for i=1:1:num_steps
        figure(fig_num);
        temp=pos_cell{i};
        scatter(temp(:,1),temp(:,2),10,[0 0 1],'o','filled');
        grid off
        xlim([-lenx lenx])
        ylim([-leny leny])
        title=(['Step =', num2str(i)]);
        xlabel('L/r0')
        ylabel('L/r0')
        frame=getframe;
        writeVideo(myObj,frame);
    end
    close(myObj);
end