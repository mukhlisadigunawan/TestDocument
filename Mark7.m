
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
nCH=3; % number of cluster head
% Number of Dead Nodes in the beggining %
dead_nodes=0;
% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=100;
sinky=300;
%%% Energy Values %%%
% Initial Energy of a Node (in Joules) %
Emin=0.5*10^-6;
Emax=1*10^-6;
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
locIDx=[];
locIDy=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creation of the Wireless Sensor Network %%%
% Plotting the WSN %
for i=1:n
    
    SN(i).id=i;	% sensor's ID number
    SN(i).x=rand(1,1)*xm;	% X-axis coordinates of sensor node
    SN(i).y=rand(1,1)*ym;	% Y-axis coordinates of sensor node
    
    locIDx(i)=SN(i).x;
    locIDy(i)=SN(i).y;
    
%     fprintf('x = %u , y = %u \n', [locIDx(i) locIDy(i)].');
    
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

% Sorting Highest Eo in the Network%





%while operating_nodes>7
 for z=1:2
    % Displays Current Round %
    rnd
    
    % Reseting Previous Amount Of Cluster Heads In the Network %
    CLheads=0;
    
    % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
    energy=0;
        
    % Cluster Heads Election %
    
    for i=1:n
        SN(i).cluster=0;    % reseting cluster in which the node belongs to
        SN(i).role=0;       % reseting node role
        SN(i).chid=0;       % reseting cluster head id
        
        if (SN(i).E>0)
            for i=1:n
                T(i)=SN(i).E;
            end
            
            
            A=sort(T,'descend'); % Creates array A containing the energy level each node,
            % sorted in an asceding order
            
            A_id(1:n)=0;
            % Creates array A_id which is sorted in a way that it's elements are
            % aligned with those of A. Contains the node ID
            % A = Energy
            % A_id = Number Id
            for i=1:n
                for j=1:n
                    if A(i)==SN(j).E
                        A_id(i)=SN(j).id;
                    end
                end
            end
            topID=[];
            topPower=[];
            for i=1:nCH
                topPower(i)=A(i);
                topID(i)=A_id(i);
                %   fprintf('ID = %u , Power = %u \n', [topID(i) topPower(i)].');
                
                
                SN(i).role=1;	% assigns the node role of acluster head
                SN(i).rn=rnd;	% Assigns the round that the cluster head was elected to the data table
                SN(i).tel=SN(i).tel + 1;
                SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2); % calculates the distance between the sink and the cluster hea
                CLheads=CLheads+1;	% sum of cluster heads that have been elected
                SN(i).cluster=CLheads; % cluster of which the node got elected to be cluster head
                CL(CLheads).x=SN(i).x; % X-axis coordinates of elected cluster head
                CL(CLheads).y=SN(i).y; % Y-axis coordinates of elected cluster head
                CL(CLheads).id=i; % Assigns the node ID of the newly elected cluster head to an array
                
                
            end
            %  fprintf('=================================\n');
            
        end
    end
    
    % Fixing the size of "CL" array %
    CL=CL(1:CLheads);
    
    
    
    % Grouping the Nodes into Clusters & caclulating the distance between node and cluster head %
    
    for i=1:n
        if  (SN(i).role==0) && (SN(i).E>0) && (CLheads>0) % if node is normal
            for m=1:CLheads
                d(m)=sqrt((CL(m).x-SN(i).x)^2 + (CL(m).y-SN(i).y)^2);
                % we calculate the distance 'd' between the sensor node that is
                % transmitting and the cluster head that is receiving with the following equation+
                % d=sqrt((x2-x1)^2 + (y2-y1)^2) where x2 and y2 the coordinates of
                % the cluster head and x1 and y1 the coordinates of the transmitting node
            end
            d=d(1:CLheads); % fixing the size of "d" array
            [M,I]=min(d(:)); % finds the minimum distance of node to CH
            [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
            SN(i).cluster=Col; % assigns node to the cluster
            SN(i).dtch= d(Col); % assigns the distance of node to CH
            SN(i).chid=CL(Col).id;
        end
    end
    
    %%=======================================================PEGASIS TEMENNYA PEGASUSNYA SAINT SAIYA BY JOCKI========================================================================
    %%%%%%%%%%%%%%%%%%%% Network Establishment Parameters %%%%%%%%%%%%%%%%%%%%
    %%% Area of Operation %%%
    % Field Dimensions in meters %
    xmPeg=200;
    ymPeg=200;
    xPeg=0; % added for better display results of the plot
    yPeg=0; % added for better display results of the plot
    % Number of Nodes in the field %
    q=0;
    % Number of Dead Nodes in the beggining %
    dead_nodesPeg=0;
    % Coordinates of the Sink (location is predetermined in this simulation) %
    sinkxPeg=100;
    sinkyPeg=300;
    %%% Energy Values %%%
    % Initial Energy of a Node (in Joules) %
    Eo=2; % units in Joules
    % Energy required to run circuity (both for transmitter and receiver) %
    EelecPeg=50*10^(-9); % units in Joules/bit
    ETxPeg=50*10^(-9); % units in Joules/bit
    ERxPeg=50*10^(-9); % units in Joules/bit
    % Transmit Amplifier Types %
    EampPeg=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
    % Data Aggregation Energy %
    EDAPeg=5*10^(-9); % units in Joules/bit
    % Size of data package %
    k=4000; % units in bits
    % Round of Operation %
    rndPeg=0;
    % Current Number of operating Nodes %
    operating_nodesPeg=nCH;
    transmissionsPeg=0;
    d(nCH,nCH)=0;
    dist_list(nCH,nCH)=0;
    temp_deadPeg=0;
    dead_nodesPeg=0;
    selectedPeg=0;
    flag1stdeadPeg=0;
    countPeg=0;
    turnPeg=0;
    temp_valPeg=0;
    cl_posPeg=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Creation of the Wireless Sensor Network %%%
    
    % Plotting the WSN %
    IDPEG=[];
    POWERPEG=[];
    res=1;  %separate LEACH and PEGASIS variable
    locTOPIDx=[];
    locTOPIDy=[];
    for i=1:n
        for j=1:nCH
            if i==topID(j)
                locTOPIDx(res)=locIDx(i);
                locTOPIDy(res)=locIDy(i);
                IDPEG(res) =i;
                POWERPEG(res)=topPower(j);
                %       fprintf('id(j) = %u , power = %u \n', [(i) topPower(j)].');
                res=res+1;
            end
        end
    end
    % for i=1:nCH
    %
    % end
    
    for i=1:n
        
        SN(i).idPeg=i;
        
        for x=1:nCH
            
            if i==IDPEG(x)
                SN(i).xPeg=locTOPIDx(x);	% X-axis coordinates of sensor node
                SN(i).yPeg=locTOPIDy(x);	% Y-axis coordinates of sensor node
                labels = {IDPEG(x)};
                labelse = {POWERPEG(x)};
                SN(x).re=IDPEG(x);
                hold on;
                
            else
                SN(i).xPeg=locIDx(i);	% X-axis coordinates of sensor node
                SN(i).yPeg=locIDy(i);	% Y-axis coordinates of sensor node
                labels = {(i)};
                labelse = {A(i)};
                hold on;
            end
        end
        
        SN(i).EPeg=Eo;     % nodes energy levels (initially set to be equal to "Eo"
        SN(i).condPeg=1;   % States the current condition of the node. when the node is operational its value is =1 and when dead =0
        SN(i).dtsPeg=0;    % nodes distance from the sink
        SN(i).rolePeg=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
        SN(i).posPeg=0;
        SN(i).closestPeg=0;
        SN(i).prevPeg=0;
        SN(i).disPeg=0;	% distance between two nodes headin towards to the cluster head from position 1
        SN(i).dis2Peg=0;   % distance between two nodes headin towards to the cluster head from position 2
        SN(i).orderPeg=0;
        SN(i).selPeg=0;    % states if the node has already operated for this round or not (if 0 then no, if 1 then yes)
        SN(i).ropPeg=0;    % number of rounds node was operational
        SN(i).telPeg=0;    % states how many times the node was elected as a Cluster Head
        %orderPeg(i)=0;
        
        if z==1
        hold on;
        figure(10)
        plot(xPeg,yPeg,xmPeg,ymPeg,SN(i).xPeg,SN(i).yPeg,'ob',sinkxPeg,sinkyPeg,'*r');
        text(SN(i).xPeg,SN(i).yPeg,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
        hold on;
        title 'Wireless Sensor Network Pegasis';
        xlabel '(m)';
        ylabel '(m)';
        end
        
        if z==2
        hold on;
        figure(11)
        plot(xPeg,yPeg,xmPeg,ymPeg,SN(i).xPeg,SN(i).yPeg,'ob',sinkxPeg,sinkyPeg,'*r');
        text(SN(i).xPeg,SN(i).yPeg,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
        hold on;
        title 'Wireless Sensor Network Pegasis';
        xlabel '(m)';
        ylabel '(m)';
        end
    end
    
    % Calculates Distance Between Each Node and the Sink (Base Station) %
    for i=1:n
        for j=1:nCH
            if  i==IDPEG(j)
                SN(i).dtsPeg=sqrt((sinkxPeg-SN(i).xPeg)^2 + (sinkyPeg-SN(i).yPeg)^2);
                SN(i).EPegsink=EelecPeg*k + EampPeg*k*(SN(i).dtsPeg)^2;
                T(i)=SN(i).dtsPeg;
            end
        end
    end
    
    
    
    APeg=sort(T,'descend'); % Creates array APeg containing the distance between each node and the sink,
    % sorted in an asceding order
    
    APeg_id(1:n)=0;
    % Creates array APeg_id which is sorted in a way that it's elements are
    % aligned with those of APeg. Contains the node ID
    for i=1:n
        for j=1:n
            if APeg(i)==SN(j).dtsPeg
                APeg_id(i)=SN(j).idPeg;
            end
        end
    end
    
    
    % Creation of d Array with shortest distances %
    
    
    for i=1:n
        SN(i).closestPeg=0;
        for j=1:nCH
            if  i==IDPEG(j)
                d(j,i)=sqrt((SN(i).xPeg-SN(j).xPeg)^2 + (SN(i).yPeg-SN(j).yPeg)^2);
                dist_list(j,i) = d(j,i);
                if d(j,i)==0
                    d(j,i)=9999;
                    dist_list(j,i) = d(j,i);
                end
            end
        end
    end
    
    
    for i=1:n
        
        [M,I]=min(d(:,i)); % finds the minimum distance of node to CH
        [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
        SN(x).closestPeg=Row; % assigns node to the cluster
        SN(x).disPeg= d(Row,i); % assigns the distance of node to CH
        
        
    end
    
    
    % Choosing furthest node from sink %
    for i=1:n
        
        if SN(APeg_id(i)).EPeg>0 && SN(APeg_id(i)).selPeg==0 && SN(APeg_id(i)).condPeg==1
            set= IDPEG(i);
            SN(set).selPeg=1;
            SN(set).posPeg=1;
            break;
            
        end
    end
    orderPeg(1)=set;
    
    tempPeg=1;
    while selectedPeg<nCH
        min_disPeg=9999;
        for i=1:n
            for x=1:nCH
                if i==IDPEG(x)
                    if  SN(i).selPeg==0
                        d=sqrt((SN(i).xPeg-SN(set).xPeg)^2 + (SN(i).yPeg-SN(set).yPeg)^2);
                        if d<min_disPeg
                            min_disPeg=d;
                            next=i;
                        end
                    end
                end
            end
        end
        selectedPeg=selectedPeg+1;
        SN(set).closestPeg=next;
        SN(set).disPeg=min_disPeg;
        SN(next).selPeg=1;
        SN(next).prevPeg=set;
        SN(next).dis2Peg=sqrt((SN(set).xPeg-SN(next).xPeg)^2 + (SN(set).yPeg-SN(next).yPeg)^2);
        plot([SN(set).xPeg SN(next).xPeg], [SN(set).yPeg SN(next).yPeg])
        hold on;
        
        set=next;
        tempPeg=tempPeg+1;
        orderPeg(tempPeg)=set;
        
    end
    
    orderPeg(nCH+1)=[];
    SN(set).posPeg=2;
    SN(set).disPeg=0;
    SN(set).closestPeg=0;
    
    for i=1:nCH
        if SN(i).closestPeg==set && SN(i).posPeg==0;
            SN(set).prevPeg=i;
            SN(set).dis2Peg=sqrt((SN(i).xPeg-SN(set).xPeg)^2 + (SN(i).yPeg-SN(set).yPeg)^2);
        end
    end
    
    energyPeg=0;
    
    for i=1:nCH
        SN(i).rolePeg=0;
    end
    
    % Cluster Head Election %
    
    mindis=sort(T,'ascend'); % Creates array APeg containing the distance between each node and the sink,
    % sorted in an asceding order
    c=1;
    mindis_id(1:n)=0;
    % Creates array APeg_id which is sorted in a way that it's elements are
    % aligned with those of APeg. Contains the node ID
    for i=1:n
        for j=1:n
            %fprintf('ID = %u /n , jarak = %u \n' , [mindis_id(i) , mindis(i)].');
            
            if mindis(i)==SN(j).dtsPeg
                mindis_id(c)=SN(j).idPeg;
                %fprintf('ID = %u \n' , [SN(j).idPeg ]');
                
                fprintf('ID = %u /n , jarak = %u \n' , [mindis_id(c) , mindis(i)].');
                c=c+1;
            end
        end
    end
    
    
    tes=1;
    cluster_head=mindis_id(tes);
    % fprintf('ID = %u /n , jarak = %u \n' , [mindis_id(i)].');
    if SN(cluster_head).condPeg==0
        while SN(cluster_head).condPeg==0
            tes=tes+1;
            cluster_head=mindis_id(tes);
        end
    end
    
    if z==1 && SN(cluster_head).condPeg==1
        SN(cluster_head).rolePeg=1;
        SN(cluster_head).telPeg=SN(cluster_head).telPeg+1;
        figure(10)
        plot(SN(cluster_head).xPeg,SN(cluster_head).yPeg,'+r')
    end
    
    if z==2 && SN(cluster_head).condPeg==1
        SN(cluster_head).rolePeg=1;
        SN(cluster_head).telPeg=SN(cluster_head).telPeg+1;
        figure(11)
        plot(SN(cluster_head).xPeg,SN(cluster_head).yPeg,'+r')
    end
    
    fprintf('Leadernode = %u \n', [mindis_id(tes)].');
    
    for i=1:nCH
        if orderPeg(i)==cluster_head
            orderPeg(i);
            cl_posPeg=i;
            break;
        end
    end
    
    for i=1:n
        for x=1:nCH
            if i==IDPEG(x)
                SN(i).E=SN(i).E-rand(1,1);
                fprintf('ID = %u \t , energy = %u \n', [i , SN(i).E].');
            end
        end
    end
    
     
    end

    %%%%%% Steady-State Phase %%%%%%
    
    
%     
%     % Energy Dissipation for normal nodes %
%     
%     for i=1:n
%         if (SN(i).cond==1) && (SN(i).role==0) && (CLheads>0)
%             if SN(i).E>0
%                 ETx= Eelec*k + Eamp * k * SN(i).dtch^2;
%                 SN(i).E=SN(i).E - ETx;
%                 energy=energy+ETx;
%                 
%                 % Dissipation for cluster head during reception
%                 if SN(SN(i).chid).E>0 && SN(SN(i).chid).cond==1 && SN(SN(i).chid).role==1
%                     ERx=(Eelec+EDA)*k;
%                     energy=energy+ERx;
%                     SN(SN(i).chid).E=SN(SN(i).chid).E - ERx;
%                     if SN(SN(i).chid).E<=0  % if cluster heads energy depletes with reception
%                         SN(SN(i).chid).cond=0;
%                         SN(SN(i).chid).rop=rnd;
%                         dead_nodes=dead_nodes +1;
%                         operating_nodes= operating_nodes - 1
%                     end
%                 end
%             end
%             
%             
%             if SN(i).E<=0       % if nodes energy depletes with transmission
%                 dead_nodes=dead_nodes +1;
%                 operating_nodes= operating_nodes - 1
%                 SN(i).cond=0;
%                 SN(i).chid=0;
%                 SN(i).rop=rnd;
%             end
%             
%         end
%     end
%     
%     
%     
%     % Energy Dissipation for cluster head nodes %
%     
%     for i=1:n
%         if (SN(i).cond==1)  && (SN(i).role==1)
%             if SN(i).E>0
%                 ETx= (Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
%                 SN(i).E=SN(i).E - ETx;
%                 energy=energy+ETx;
%             end
%             if  SN(i).E<=0     % if cluster heads energy depletes with transmission
%                 dead_nodes=dead_nodes +1;
%                 operating_nodes= operating_nodes - 1
%                 SN(i).cond=0;
%                 SN(i).rop=rnd;
%             end
%         end
%     end
%     
%     
%     if operating_nodes<n && temp_val==0
%         temp_val=1;
%         flag1stdead=rnd
%     end
%     % Display Number of Cluster Heads of this round %
%     %CLheads;
%     
%     
%     transmissions=transmissions+1;
%     if CLheads==0
%         transmissions=transmissions-1;
%     end
%     
%     
%     % Next Round %
%     rnd= rnd +1;
%     
%     tr(transmissions)=operating_nodes;
%     op(rnd)=operating_nodes;
%     
%     if energy>0
%         nrg(transmissions)=energy;
%     end
%     
%     
% end
% sum=0;
% for i=1:flag1stdead
%     sum=nrg(i) + sum;
% end
% temp1=sum/flag1stdead;
% temp2=temp1/n;
% for i=1:flag1stdead
%     avg_node(i)=temp2;
% end

% % Plotting Simulation Results "Operating Nodes per Round" %
% figure(2)
% plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
% title ({'LEACH'; 'Operating Nodes per Round';})
% xlabel 'Rounds';
% ylabel 'Operational Nodes';
% hold on;
% 
% % Plotting Simulation Results  %
% figure(3)
% plot(1:transmissions,tr(1:transmissions),'-r','Linewidth',2);
% title ({'LEACH'; 'Operational Nodes per Transmission';})
% xlabel 'Transmissions';
% ylabel 'Operational Nodes';
% hold on;
% 
% % Plotting Simulation Results  %
% figure(4)
% plot(1:flag1stdead,nrg(1:flag1stdead),'-r','Linewidth',2);
% title ({'LEACH'; 'Energy consumed per Transmission';})
% xlabel 'Transmission';
% ylabel 'Energy ( J )';
% hold on;
% 

