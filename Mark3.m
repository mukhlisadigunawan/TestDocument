%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          LEACH Protocol                              %
%                          Part of Thesis:                             %
%       Energy-Efficient Protocols In Wireless Sensor Networks         %
%                                                                      %
% (c) Alexandros Nikolaos Zattas                                       %
% University of Peloponnese                                            %
% Department of Informatics and Telecommunications                     %
% For any related questions or additional information                  %
% send an e-mail to:                                                   %
% alexzattas@gmail.com                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;    
%%%%%%%%%%%%%%%%%%%% Network Establishment Parameters %%%%%%%%%%%%%%%%%%%%
%%% Area of Operation %%%
% Field Dimensions in meters %
xm=200;
ym=200;
x=0; % added for better display results of the plot
y=0; % added for better display results of the plot
% Number of Nodes in the field %
n=10;
% Number of Dead Nodes in the beggining %
dead_nodes=0;
% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=100;
sinky=100;
%%% Energy Values %%%
% Initial Energy of a Node (in Joules) % 
Emin=0.5;
Emax=1;
%Eo=(Eob-Eoa).*rand(1,1)+Eoa; % units in Joules
% Energy required to run circuity (both for transmitter and receiver) %
Eelec=50*10^(-9); % units in Joules/bit
%ETx=50*10^(-9); % units in Joules/bit
%ERx=10*10^(-9); % units in Joules/bit
% Transmit Amplifier Types %
Eamp=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
% Data Aggregation Energy %
EDA=5*10^(-9); % energy for data aggregation(Joules/bit)
% Size of data package %
k=2000; % units in bits
% Suggested percentage of cluster head %
p=0.05;% a 5 percent of the total amount of nodes used in the network is proposed to give good results
% Number of Clusters %
No=p*n; 
% Round of Operation %
rnd=0;
% Current Number of operating Nodes %
operating_nodes=n;
transmissions=0;
temp_val=0;
flag1stdead=0;

%PEGASIS%
d(n,n)=0;
temp_dead=0; 
dead_nodes=0;
selected=0; 
count=0; 
turn=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Creation of the Wireless Sensor Network %%%
% Plotting the WSN %
for i=1:n
    
    SN(i).id=i;	% sensor's ID number
    SN(i).x=rand(1,1)*xm;	% X-axis coordinates of sensor node
    SN(i).y=rand(1,1)*ym;	% Y-axis coordinates of sensor node
    SN(i).E=(Emax-Emin).*rand(1,1)+Emin;     % nodes energy levels (initially set to be equal to "Eo"
    SN(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    SN(i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
    SN(i).rop=0;	% number of rounds node was operational
    SN(i).dts=0;    % nodes distance from the sink
    SN(i).tel=0;	% states how many times the node was elected as a Cluster Head
    
    %LEACH%
    SN(i).rn=0;     % round node got elected as cluster head
    SN(i).chid=0;   % node ID of the cluster head which the "i" normal node belongs to
    SN(i).cluster=0;	% the cluster which a node belongs to
    SN(i).rleft=0;  % rounds left for node to become available for Cluster Head election
    SN(i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs

    %PEGASIS%
    SN(i).pos=0;
    SN(i).closest=0;
    SN(i).prev=0;
    SN(i).dis=0;	% distance between two nodes headin towards to the cluster head from position 1
    SN(i).dis2=0;   % distance between two nodes headin towards to the cluster head from position 2
    SN(i).order=0;
    SN(i).sel=0;    % states if the node has already operated for this round or not (if 0 then no, if 1 then yes) 
    order(i)=0;
    
    labels = {SN(i).id};
    hold on;
    figure(1)
    plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r');
    text(SN(i).x,SN(i).y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    
end
 
                      %%%%%% Set-Up Phase %%%%%% 
                      
             
while operating_nodes>0
        
    % Displays Current Round %     
    rnd     
	% Threshold Value %
	t=(p/(1-p*(mod(rnd,1/p))));
    
    % Re-election Value %
    tleft=mod(rnd,1/p);
 
	% Reseting Previous Amount Of Cluster Heads In the Network %
	CLheads=0;
    
    % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
    energy=0;
 
    
          
    % Cluster Heads Election %
         for i=1:n
            SN(i).cluster=0;    % reseting cluster in which the node belongs to
            SN(i).role=0;       % reseting node role
            SN(i).chid=0;       % reseting cluster head id
            if SN(i).rleft>0
               SN(i).rleft=SN(i).rleft-1;
            end
            if (SN(i).E>0) && (SN(i).rleft==0)
                generate=rand;	
                    if generate< t
                    SN(i).role=1;	% assigns the node role of acluster head
                    SN(i).rn=rnd;	% Assigns the round that the cluster head was elected to the data table
                    SN(i).tel=SN(i).tel + 1;   
                    SN(i).rleft=1/p-tleft;    % rounds for which the node will be unable to become a CH
                    SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2); % calculates the distance between the sink and the cluster hea
                    CLheads=CLheads+1;	% sum of cluster heads that have been elected 
                    SN(i).cluster=CLheads; % cluster of which the node got elected to be cluster head
                    CL(CLheads).x=SN(i).x; % X-axis coordinates of elected cluster head
                    CL(CLheads).y=SN(i).y; % Y-axis coordinates of elected cluster head
                    CL(CLheads).id=i; % Assigns the node ID of the newly elected cluster head to an array
                    end
        
            end
        end
        
	% Fixing the size of "CL" array %
	CL=CL(1:CLheads);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PEGASIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   % Calculates Distance Between Each Node and the Sink (Base Station) %
 for i=1:n
    SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2);
    SN(i).Esink=Eelec*k + Eamp*k*(SN(i).dts)^2;
    T(i)=SN(i).dts;
 end
 
 
    A=sort(T,'descend'); % Creates array A containing the distance between each node and the sink,
                % sorted in an asceding order
     
     A_id(1:n)=0;
     % Creates array A_id which is sorted in a way that it's elements are
     % aligned with those of A. Contains the node ID
     for i=1:n
         for j=1:n
            if A(i)==SN(j).dts
               A_id(i)=SN(j).id;
            end
         end
     end
     
     
     % Creation of d Array with shortest distances %
     
     
            for i=1:n
             SN(i).closest=0;
             for j=1:n
                d(j,i)=sqrt((SN(i).x-SN(j).x)^2 + (SN(i).y-SN(j).y)^2);
                if d(j,i)==0
                    d(j,i)=9999;
                end
             end
            end
       
                  
        for i=1:n     
            [M,I]=min(d(:,i)); % finds the minimum distance of node to CH
            [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
            SN(i).closest=Row; % assigns node to the cluster
            SN(i).dis= d(Row,i); % assigns the distance of node to CH
        end
     
     
        % Choosing furthest node from sink %
        for i=1:n
             if SN(A_id(i)).E>0 && SN(A_id(i)).sel==0 && SN(A_id(i)).cond==1
                set= A_id(i);
                SN(set).sel=1;
                SN(set).pos=1;
                break;
             end
        end
     order(1)=set;

     temp=1;   
        while selected<n
            min_dis=9999;
            for i=1:n
                if  SN(i).sel==0 
                    d=sqrt((SN(i).x-SN(set).x)^2 + (SN(i).y-SN(set).y)^2);
                    if d<min_dis
                        min_dis=d;
                        next=i; 
                    end
                end
            end
            selected=selected+1;
            SN(set).closest=next;
            SN(set).dis=min_dis;
            SN(next).sel=1;
            SN(next).prev=set;
            SN(next).dis2=sqrt((SN(set).x-SN(next).x)^2 + (SN(set).y-SN(next).y)^2);
            plot([SN(set).x SN(next).x], [SN(set).y SN(next).y])
            hold on;
            set=next;
            temp=temp+1;
            order(temp)=set;
        end
        
   
        order(n+1)=[];
        SN(set).pos=2;
        SN(set).dis=0;
        SN(set).closest=0;
        for i=1:n
            if SN(i).closest==set && SN(i).pos==0;
               SN(set).prev=i;
               SN(set).dis2=sqrt((SN(i).x-SN(set).x)^2 + (SN(i).y-SN(set).y)^2);
            end
        end
       
         
       
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEACH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                           %%%%%% Steady-State Phase %%%%%%
                      
  
% Energy Dissipation for normal nodes %
    
    for i=1:n
       if (SN(i).cond==1) && (SN(i).role==0) && (CLheads>0)
       	if SN(i).E>0
            ETx= Eelec*k + Eamp * k * SN(i).dtch^2;
            SN(i).E=SN(i).E - ETx;
            energy=energy+ETx;
            
        % Dissipation for cluster head during reception
        if SN(SN(i).chid).E>0 && SN(SN(i).chid).cond==1 && SN(SN(i).chid).role==1
            ERx=(Eelec+EDA)*k;
            energy=energy+ERx;
            SN(SN(i).chid).E=SN(SN(i).chid).E - ERx;
             if SN(SN(i).chid).E<=0  % if cluster heads energy depletes with reception
                SN(SN(i).chid).cond=0;
                SN(SN(i).chid).rop=rnd;
                dead_nodes=dead_nodes +1;
                operating_nodes= operating_nodes - 1
             end
        end
        end
        
        
        if SN(i).E<=0       % if nodes energy depletes with transmission
        dead_nodes=dead_nodes +1;
        operating_nodes= operating_nodes - 1
        SN(i).cond=0;
        SN(i).chid=0;
        SN(i).rop=rnd;
        end
        
       end
        end            
    
    
    
% Energy Dissipation for cluster head nodes %
   
   for i=1:n
     if (SN(i).cond==1)  && (SN(i).role==1)
         if SN(i).E>0
            ETx= (Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
            SN(i).E=SN(i).E - ETx;
            energy=energy+ETx;
         end
         if  SN(i).E<=0     % if cluster heads energy depletes with transmission
         dead_nodes=dead_nodes +1;
         operating_nodes= operating_nodes - 1
         SN(i).cond=0;
         SN(i).rop=rnd;
         end
     end
   end
   
  
    if operating_nodes<n && temp_val==0
        temp_val=1;
        flag1stdead=rnd
    end
    % Display Number of Cluster Heads of this round %
    %CLheads;
   
    
    transmissions=transmissions+1;
    if CLheads==0
    transmissions=transmissions-1;
    end
    
 
    % Next Round %
    rnd= rnd +1;
    
    tr(transmissions)=operating_nodes;
    op(rnd)=operating_nodes;
    
    if energy>0
    nrg(transmissions)=energy;
    end
    
end
sum=0;
for i=1:flag1stdead
    sum=nrg(i) + sum;
end
temp1=sum/flag1stdead;
temp2=temp1/n;
for i=1:flag1stdead
avg_node(i)=temp2;
end
    
    % Plotting Simulation Results "Operating Nodes per Round" %
    figure(2)
    plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
    title ({'LEACH'; 'Operating Nodes per Round';})
    xlabel 'Rounds';
    ylabel 'Operational Nodes';
    hold on;
    
    % Plotting Simulation Results  %
    figure(3)
    plot(1:transmissions,tr(1:transmissions),'-r','Linewidth',2);
    title ({'LEACH'; 'Operational Nodes per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Operational Nodes';
    hold on;
    
    % Plotting Simulation Results  %
    figure(4)
    plot(1:flag1stdead,nrg(1:flag1stdead),'-r','Linewidth',2);
    title ({'LEACH'; 'Energy consumed per Transmission';})
    xlabel 'Transmission';
    ylabel 'Energy ( J )';
    hold on;    
    
    % Plotting Simulation Results  %
    %figure(5)
    %plot(1:flag1stdead,avg_node(1:flag1stdead),'-r','Linewidth',2);
    %title ({'LEACH'; 'Average Energy consumed by a Node per Transmission';})
    %xlabel 'Transmissions';
    %ylabel 'Energy ( J )';
    %hold on;