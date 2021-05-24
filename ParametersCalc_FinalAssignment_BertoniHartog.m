% Parameters
N = 1000; % Number of mussels
Length = 50; % Length of the Arena, in cm
EndTime = 500; % Number of minutes
dT = 1; % minute
DeltaX = 1;
DeltaY = 1;
Nsim = 1;

% Declaring the parameters that describe movement speed
P1 = 100; 

P3 = 3; 
D1 = 1; % cm, mussels within this distance are considered neighbors
D2 = 4; % cm, mussels within this distance are considered as "within cluster distance", in cm.
 
% Parameters related to phytoplankton
% massPhyto = 2.4; % microg
% gamma = 0.1;  % per minute 0.3/60/24
% k = 108; %micro g per minute 108*10^(-6)
D = 2.1*10^(-5)*10000;     % cm^2 per minute 0.03/60/24
F = 24.26*10^3; % (cell/ml) average phytoplankton eaten by 1mussel per timestep
AveragePhytoDensity = 2256*10^3;     % average phytoplankton density in sea water [cell/ml]   
Rescale = -80*AveragePhytoDensity; % This parameter will multiply the phytoplankton concentration (P). It is inserted so that the step size resemble on average the one used in the old model



D3 = 2;
D4 = 6;
density = 0.5;
radius = linspace(0,10,50);
max_mussel_within_density = zeros(1,50);

Cluster3 = zeros(EndTime,Nsim);
Cluster4 = zeros(EndTime, Nsim);




%% Various simulation loop
for sim=1:Nsim
    % Constant and variable declaration, used below
    Diagonal=diag(ones(N,1)); % Used below
    Distance=zeros(N,N); % Assigning the array
    Beta = zeros(N,1);
    StepSize = zeros(N,1);
    Angle = zeros(N,1);


    % Initial distribution
    X = rand(N,1)*50;
    Y = rand(N,1)*50;
    X_int = X;
    Y_int = Y;
    P = ones(Length, Length)*AveragePhytoDensity;  % Phytoplankton distribution
    dP = zeros(Length, Length);  % Change in phytoplankton concentration
    NetP = zeros(Length, Length);   % Phytoplankton net flux of diffusion 
    
    %% Main loop
  
    for Time=1:EndTime

        DistributionMussels =zeros(Length,Length);

        X_int = round(X);
        Y_int = round(Y);

       % Finding the amount of mussels in each cell
       for m = 1:Length
           for n= 1:Length
               for s=1:N
                   if X_int(s)+1==m && Y_int(s)+1==n
                       DistributionMussels(m,n) = DistributionMussels(m,n) + 1;
                   end 
               end
           end
       end



       % calculate the reaction term
       dP = -DistributionMussels*F;


       % Calculate the diffusion flows
       FXP (1:Length,2:Length) = -D * (P(1:Length,2:Length)- P(1:Length,1:Length-1))/DeltaX;
       FYP (2:Length,1:Length) = -D * (P(2:Length,1:Length)- P(1:Length-1,1:Length))/DeltaY;

       FXP(1:Length,1) = -D * (P(1:Length,1)- P(1:Length,Length))/DeltaX;       % Boundaries conditions for the flows
       FYP(1, 1:Length) = -D * (P(1,1:Length)- P(Length,1:Length))/DeltaY;

       for i=1:Length
           for j=1:Length
               % boundaries conditions
                if i~=Length && j==Length  
                   NetP(i,j) = FXP(i,j)-FXP(i,1)+FYP(i, j) - FYP(i+1,j);
                elseif i==Length && j==Length
                   NetP(i,j) = FXP(i,j)-FXP(i,1)+FYP(i, j) - FYP(1,j);
                elseif i==Length && j~=Length  
                   NetP(i, j) = FXP(i, j)-FXP(i,j+1)+FYP(i, j) - FYP(1,j);
                % regular flow
                else 
                   NetP(i,j)= FXP(i,j)-FXP(i,j+1)+ FYP(i,j)-FYP(i+1,j);
                end

           end
       end

       % Updating the phytoplankton population

       P = P+NetP*dT+dP*dT.*P/AveragePhytoDensity; % *


       % MUSSEL MOVEMENT CALCULATIONS

       for i=1:N
         for j=1:N
         Distance(i,j)= sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2);
         end
       end

        Nr_In_Dist1 = (Distance<D1)-Diagonal;
        
        Nr_In_Dist3 = (Distance<D3)-Diagonal;
        Nr_In_Dist4 = (Distance<D4)-Diagonal;
        
        % Calculating densities of mussles
        V1 = sum(Nr_In_Dist1)/(D1.^2*pi()); % sum/surface
        
        V3 = sum(Nr_In_Dist3)/(D3.^2*pi()); % sum/surface
        V4 = sum(Nr_In_Dist4)/(D4.^2*pi()); % sum/surface
        
        
        Cluster3(Time,sim) = sum(V3>density);
        Cluster4(Time,sim) = sum(V4>density);
        

        
        
        for m = 1:Length
           for n= 1:Length
               for s=1:N
                   if X_int(s)+1==m && Y_int(s)+1==n
                       Beta(s)=1/(max(0.001,P1*V1(s)-0.34*(80+10*mean(mean(P))/P(m,n)))+P3);
                       % Setting a random step size
                       StepSize(s)= -Beta(s)* log(rand());

                   end 
               end
           end
        end

        % setting a random stepsize
        Angle = rand(N,1)*360;

        %Implementing the steps
        X = X + StepSize.*cos(Angle);
        Y = Y + StepSize.*sin(Angle);

        % Periodic boundary conditions
        for i = 1:N
            if X(i)>Length
                X(i) = X(i) - Length;
            end

            if X(i)<0
                X(i) = X(i) + Length;
            end

            if Y(i)>Length
                Y(i) = Y(i) - Length;
            end

            if Y(i)<0
                Y(i) = Y(i) + Length;
            end
        end

%         % Mussel distribution
%         figure(1)
%         plot(X,Y,'g.','MarkerSize',20)
%         title('Mussels distribution, N=1000')
%         ylabel('Y (cm)')
%         xlabel('X (cm)')
%         ylim([0 50])
%         xlim([0 50])
%         drawnow

%         % food distribution
%         figure(2) 
%         imagesc(P)
%         colorbar()
%         title('Phytoplankton distribution')
%         ylabel('Y (cm)')
%         xlabel('X (cm)')
%         ylim([0 50])
%         xlim([0 50])

    %     % Mussel and food distribution 
    %     figure(3) 
    %     imagesc(P)
    %     colorbar()
    %     hold on
    %     plot(X,Y,'g.','MarkerSize',20)
    %     title('Mussels and food distribution, N=1000')
    %     ylabel('Y (cm)')
    %     xlabel('X (cm)')
    %     ylim([0 50])
    %     xlim([0 50])
    %     drawnow



    end
end  
 
%%

for i=1:50
    Nr_In_Dist = (Distance<radius(i))-Diagonal;
    V = sum(Nr_In_Dist)/(radius(i).^2*pi());
    Cluster = sum(V>density);
    max_mussel_within_density(i) = max(Cluster);
end 

figure()
plot(radius, max_mussel_within_density/N*100,'g-','LineWidth', 2)
hold on 
plot(radius,old_max_mussel_within_density/N*100,'b-','LineWidth', 2)
xlabel('Radius [cm]','Interpreter','latex','FontSize',14)
ylabel('\% of mussels','Interpreter','latex','FontSize',14)
title('Mussel with a local density of 0.5 $[muss/cm^2]$ ','Interpreter','latex','FontSize',12)
legend('With food', 'No food','Interpreter','latex','FontSize',12)

figure()
for p=1:Nsim
plot(Cluster3(:,p)/N*100,'g+', 'MarkerSize', 5)
hold on
plot(Cluster4(:,p)/N*100,'g-', 'MarkerSize', 4)
hold on
plot(old_Cluster3(:,p)/N*100,'b+', 'MarkerSize', 5)
hold on
plot(old_Cluster4(:,p)/N*100,'b-', 'MarkerSize', 4)
xlabel('Time [min]','Interpreter','latex','FontSize',14)
ylabel('\% of mussels','Interpreter','latex','FontSize',14)
title('Mussel aggregation over time ','Interpreter','latex','FontSize',12)
legend('With food, D=1 cm','With food, D=4 cm', 'No food, D=1 cm', 'No food, D= 4cm','Interpreter','latex','FontSize',10)
end
%     
    
    
    
    