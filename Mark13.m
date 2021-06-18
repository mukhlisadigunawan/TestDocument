% inisialisasi node
close all;
clear;
clc;
xm=200;             %ukuran dimensi plot
ym=200;             %ukuran dimensi plot
x=0;
y=0;
n=6;               %jumlah node yang di inisialisasi
p=0.05;             %faktor pengali untuk menentukan cluster head
No=p*n;             %jumlah cluster head dengan faktor pengali P
% nCH=5;              % jumlah CH yang dijadikan CH pegasis
dead_nodes=0;       %jumlah node yang mati
sinkx=100;          %koordinat Sink Node
sinky=300;          %koordinat Sink Node
Emin=0.5;           %inisialisai energi minum
Emax=1;             %inisialisasi energi maksimum
Eelec=50*10^(-9);   %energi yang dibutuhkan untuk running ( transmiter atau reciver)
Eamp=100*10^(-12);  %faktor pengali amplifire
EDA=5*10^(-9);      %energi agregasi data
k=4000;             %bersar paket yangdikirimkan
rnd=0;              %ronde dari operasi
operating_nodes=n;  %jumlah operasi node
transmissions=0;    %variable transmisi
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

DN=[];
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
    SN(i).posPeg=0;
    SN(i).closest=0;
    SN(i).prev=0;
    SN(i).dis=0;	% distance between two nodes headin towards to the cluster head from position 1
    SN(i).dis2=0;   % distance between two nodes headin towards to the cluster head from position 2
    SN(i).orderPeg=0;
    SN(i).sel=0;    % states if the node has already operated for this round or not (if 0 then no, if 1 then yes)
    orderPeg(i)=0;
    
    labels = {SN(i).id};
    hold on;
    figure(1)
    plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r');
    text(SN(i).x,SN(i).y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
    title 'WSN LEACH';
    xlabel '(m)';
    ylabel '(m)';
    figure(2)
    plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r');
    text(SN(i).x,SN(i).y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
    title 'WSN LEACH PEGASIS';
    xlabel '(m)';
    ylabel '(m)';
    
end

%%%%%% Set-Up Phase %%%%%%

while operating_nodes>0
    % for z=1:10
    nCH=ceil(operating_nodes*0.2);
    rnd;
    CLheads=0;
    energy=0;
    % Cluster Heads Election %
    for i=1:n
        SN(i).cluster=0;    % reseting cluster in which the node belongs to
        SN(i).role=0;       % reseting node role
        SN(i).chid=0;       % reseting cluster head id
    end
    %     if (SN(i).E>0)
    for i=1:n
        T(i)=SN(i).E;
    end
    A=sort(T,'descend'); % Creates array A containing the energy level each node,
    A_id(1:n)=0;
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
        SN(topID(i)).role=1;	% assigns the node role of acluster head
        SN(topID(i)).rn=rnd;	% Assigns the round that the cluster head was elected to the data table
        SN(topID(i)).tel=SN(topID(i)).tel + 1;
        SN(topID(i)).dts=sqrt((sinkx-SN(topID(i)).x)^2 + (sinky-SN(topID(i)).y)^2); % calculates the distance between the sink and the cluster hea
        CLheads=CLheads+1;	% sum of cluster heads that have been elected
        SN(topID(i)).cluster=CLheads; % cluster of which the node got elected to be cluster head
        CL(CLheads).x=SN(topID(i)).x; % X-axis coordinates of elected cluster head
        CL(CLheads).y=SN(topID(i)).y; % Y-axis coordinates of elected cluster head
        CL(CLheads).id=topID(i); % Assigns the node ID of the newly elected cluster head to an array
    end
    %     end
    % Fixing the size of "CL" array %
    CL=CL(1:CLheads);
    % Grouping the Nodes into Clusters & caclulating the distance between node and cluster head %
    for i=1:n
        %         fprintf('ID = %u , Role = %u \n', [(i) SN(i).role].');
        if  (SN(i).role==0) && (SN(i).E>0) && (CLheads>0) % if node is normal
            for m=1:CLheads
                d(m)=sqrt((CL(m).x-SN(i).x)^2 + (CL(m).y-SN(i).y)^2);
            end
            d=d(1:CLheads); % fixing the size of "d" array
            [M,I]=min(d(:)); % finds the minimum distance of node to CH
            [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
            SN(i).cluster=Col; % assigns node to the cluster
            SN(i).dtch= d(Col); % assigns the distance of node to CH
            SN(i).chid=CL(Col).id;
            if rnd==1
                figure(1)
                 plot([SN(i).x SN(SN(i).chid).x], [SN(i).y SN(SN(i).chid).y])
                figure(2)
            plot([SN(i).x SN(SN(i).chid).x], [SN(i).y SN(SN(i).chid).y])
            end
                        fprintf('ID = %u \t , ch = %u \n', [i , SN(i).chid].');
        end
    end
    
    
    q=0;
    dead_nodesPeg=0;
    EelecPeg=50*10^(-9); % units in Joules/bit
    ETxPeg=50*10^(-9); % units in Joules/bit
    ERxPeg=50*10^(-9); % units in Joules/bit
    EampPeg=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
    EDAPeg=5*10^(-9); % units in Joules/bit
    rndPeg=0;
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
                
                SN(i).EPeg=SN(i).E;     % nodes energy levels (initially set to be equal to "Eo"
                SN(i).condPeg=SN(i).cond;   % States the current condition of the node. when the node is operational its value is =1 and when dead =0
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
                hold on;
            else
                SN(i).xPeg=locIDx(i);	% X-axis coordinates of sensor node
                SN(i).yPeg=locIDy(i);	% Y-axis coordinates of sensor node
                labels = {(i)};
                labelse = {A(i)};
                hold on;
            end
        end
        
        
        orderPeg(i)=0;
        
        if rnd==1
            hold on;
            figure(3)
            plot(x,y,xm,ym,SN(i).xPeg,SN(i).yPeg,'ob',sinkx,sinky,'*r');
            text(SN(i).xPeg,SN(i).yPeg,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
            hold on;
            title 'WSN PEGASIS';
            xlabel '(m)';
            ylabel '(m)';
        end
%         
%         if rnd==1
%             hold on;
%             figure(11)
%             plot(x,y,xm,ym,SN(i).xPeg,SN(i).yPeg,'ob',sinkx,sinky,'*r');
%             text(SN(i).xPeg,SN(i).yPeg,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
%             hold on;
%             title 'Wireless Sensor Network Pegasis';
%             xlabel '(m)';
%             ylabel '(m)';
%         end
    end
    
    % Calculates Distance Between Each Node and the Sink (Base Station) %
    for i=1:n
        for j=1:nCH
            if  i==IDPEG(j)
                SN(i).dtsPeg=sqrt((sinkx-SN(i).xPeg)^2 + (sinky-SN(i).yPeg)^2);
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
    
    
    for i=1:nCH
        SN(i).closestPeg=0;
        for j=1:nCH
            %             for x=1:nCH
            %                 if  i==IDPEG(x)
            d(j,i)=sqrt((SN(i).xPeg-SN(j).xPeg)^2 + (SN(i).yPeg-SN(j).yPeg)^2);
            dist_list(j,i) = d(j,i);
            if d(j,i)==0
                d(j,i)=9999;
                dist_list(j,i) = d(j,i); %Creates array containing distance between nodes
                %                     end
                %                 end
            end
        end
    end
    
    
    for i=1:nCH
        
        [M,I]=min(d(:,i)); % finds the minimum distance of node to CH
        %         fprintf('size d = %u, ',[size(d)].');
        [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
        SN(IDPEG(i)).closestPeg=IDPEG(Row); % assigns node to the cluster
        %         fprintf(' d = %u, ',[d(Row,i)].');
        SN(IDPEG(i)).disPeg= d(Row,i); % assigns the distance of node to CH
        %         fprintf('id = %u , closest = %u, dis = %u \n' , [(IDPEG(i)) SN(IDPEG(i)).closestPeg , SN(IDPEG(i)).disPeg].');
        
        
    end
    
    
    % Choosing furthest node from sink %
    for i=1:nCH
        
        if SN(APeg_id(i)).EPeg>0 && SN(APeg_id(i)).selPeg==0 && SN(APeg_id(i)).condPeg==1
            set= APeg_id(i);
            SN(set).selPeg=1;
            SN(set).posPeg=1;
            %             fprintf('id = %u' , [ set].');
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
        if rnd==1
            figure(2)
            plot([SN(set).xPeg SN(next).xPeg], [SN(set).yPeg SN(next).yPeg])
            figure(3)
            plot([SN(set).xPeg SN(next).xPeg], [SN(set).yPeg SN(next).yPeg])
        end
        
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
        if SN(i).closestPeg==set && SN(i).posPeg==0
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
                
                %                 fprintf('ID = %u /n , jarak = %u \n' , [mindis_id(c) , mindis(i)].');
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
    
    if  SN(cluster_head).condPeg==1
        SN(cluster_head).rolePeg=1;
        SN(cluster_head).telPeg=SN(cluster_head).telPeg+1;
        if rnd==1
            plot(SN(cluster_head).xPeg,SN(cluster_head).yPeg,'+r')
        end
    end
    
    %     if z==2 && SN(cluster_head).condPeg==1
    %         SN(cluster_head).rolePeg=1;
    %         SN(cluster_head).telPeg=SN(cluster_head).telPeg+1;
    %         figure(11)
    %         plot(SN(cluster_head).xPeg,SN(cluster_head).yPeg,'+r')
    %     end
    
    %     fprintf('Leadernode = %u \n', [mindis_id(tes)].');
    
    
    
    %%%%%% Steady-State Phase %%%%%%
    %      fprintf('------------------------ \n');
    
    %LEACH Transmission%
    for i=1:n
        if (SN(i).cond==1) && (SN(i).role==0) && (CLheads>0)
            % Member Transmit %
            if SN(i).E>0
                ETx= Eelec*k + Eamp * k * SN(i).dtch^2;
                SN(i).E=SN(i).E - ETx;
                energy=energy+ETx;
                
                %                 fprintf('node ID transmit LEACH= %u \t , energy = %u \n', [i , SN(i).E].');
                
                % CH Receive%
                if SN(SN(i).chid).E>0 && SN(SN(i).chid).cond==1 && SN(SN(i).chid).role==1
                    ERx=(Eelec+EDA)*k;
                    energy=energy+ERx;
                    SN(SN(i).chid).E=SN(SN(i).chid).E - ERx;
                    
                    if SN(SN(i).chid).E<=0  % if cluster heads energy depletes with reception
                        SN(SN(i).chid).cond=0;
                        SN(SN(i).chid).rop=rnd;
                        dead_nodes=dead_nodes +1;
                        operating_nodes= operating_nodes - 1;
                    end
                    %                     fprintf('CH receive LEACH= %u \t , energy = %u \n', [SN(i).chid , SN(SN(i).chid).E].');
                    
                    
                end
                
                
            end
            %                          fprintf('------------------------ \n');
            
        end
    end
    
    %     fprintf('---------PEGASIS--------------- \n');
    
    % PEGASIS Transmission %
    for i=1:nCH
        if orderPeg(i)==cluster_head
            %             orderPeg(i);
            cl_pos=i;
            %                                     fprintf('Leadernode = %u \n', [cluster_head].');
            %             fprintf('CH = %u \n', [cl_pos].');
            break;
        end
    end
    
    for i=1:nCH
        %         fprintf('CH= %u \t , role = %u \n', [SN(orderPeg(i)).id , SN(orderPeg(i)).role].');
        
        if SN(orderPeg(i)).E>0 && SN(orderPeg(i)).cond==1  && (SN(orderPeg(i)).role==1)
            if i<cl_pos
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if SN(orderPeg(i)).posPeg==1 && SN(orderPeg(i)).rolePeg==0
                    ETx= Eelec*k + Eamp*k*SN(orderPeg(i)).dis^2;
                    SN(orderPeg(i)).E=SN(orderPeg(i)).E-ETx;
                    energy=energy+ETx;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if SN(orderPeg(i)).posPeg==0 && SN(orderPeg(i)).rolePeg==0
                    ERx=(EDA+Eelec)*k;
                    ETx= (EDA+Eelec)*k + Eamp*k*SN(orderPeg(i)).dis^2;
                    SN(orderPeg(i)).E=SN(orderPeg(i)).E-ETx-ERx;
                    energy=energy+ETx+ERx;
                end
                %                                 fprintf('ID = %u \t , energy = %u \n', [SN(orderPeg(i)).id , SN(orderPeg(i)).E].');
                %                                 fprintf('dis1 id= %u \n', [SN(orderPeg(i)).id].');
            end
            
            if i>cl_pos
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if SN(orderPeg(i)).posPeg==2 && SN(orderPeg(i)).rolePeg==0
                    ETx= Eelec*k + Eamp*k*SN(orderPeg(i)).dis2^2;
                    SN(orderPeg(i)).E=SN(orderPeg(i)).E-ETx;
                    energy=energy+ETx;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if SN(orderPeg(i)).posPeg==0 && SN(orderPeg(i)).rolePeg==0
                    ERx=(EDA+Eelec)*k;
                    ETx= (EDA+Eelec)*k + Eamp*k*SN(orderPeg(i)).dis2^2;
                    SN(orderPeg(i)).E=SN(orderPeg(i)).E-ETx-ERx;
                    energy=energy+ETx+ERx;
                end
                %                                 fprintf('ID = %u \t , energy = %u \n', [SN(orderPeg(i)).id , SN(orderPeg(i)).E].');
                %                                 fprintf('dis1 id= %u \n', [SN(orderPeg(i)).id].');
            end
            
            %Leader Node%
            if i==cl_pos
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if SN(orderPeg(i)).posPeg==0  && SN(orderPeg(i)).rolePeg==1
                    ERx=(EDA+Eelec)*k*2;
                    ETx= (EDA+Eelec)*k + Eamp*k*SN(orderPeg(i)).dts^2;
                    SN(orderPeg(i)).E=SN(orderPeg(i)).E-ETx-ERx;
                    energy=energy+ETx+ERx;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if SN(orderPeg(i)).posPeg==1  && SN(orderPeg(i)).rolePeg==1
                    ERx=(EDA+Eelec)*k;
                    ETx= (EDA+Eelec)*k + Eamp*k*SN(orderPeg(i)).dts^2;
                    SN(orderPeg(i)).E=SN(orderPeg(i)).E-ETx-ERx;
                    energy=energy+ETx+ERx;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if  SN(orderPeg(i)).posPeg==2  && SN(orderPeg(i)).rolePeg==1
                    ERx=(EDA+Eelec)*k;
                    ETx= (EDA+Eelec)*k + Eamp*k*SN(orderPeg(i)).dts^2;
                    SN(orderPeg(i)).E=SN(orderPeg(i)).E-ETx-ERx;
                    energy=energy+ETx+ERx;
                end
                %                                 fprintf('ID = %u \t , energy = %u \n', [SN(orderPeg(i)).id , SN(orderPeg(i)).E].');
                %                                 fprintf('LN id= %u \n', [SN(orderPeg(i)).id].');
            end
        end
        
        % DEAD NODES%
        if SN(orderPeg(i)).E<=0 && SN(orderPeg(i)).condPeg==1
            SN(orderPeg(i)).condPeg=0;
            SN(orderPeg(i)).cond=0;
            operating_nodes=operating_nodes-1;
            dead_nodes=dead_nodes+1;
            SN(orderPeg(i)).closestPeg=0;
            SN(orderPeg(i)).prevPeg=0;
            SN(orderPeg(i)).disPeg=0;
            SN(orderPeg(i)).dis2Peg=0;
            SN(orderPeg(i)).posPeg=101;
            SN(orderPeg(i)).ropPeg=SN(orderPeg(i)).rop;
            SN(orderPeg(i)).rop=rnd;
            %             fprintf('ID = %u \t , energy = %u \n', [SN(orderPeg(i)).id , SN(orderPeg(i)).E].');
        end
    end
    
    for i=1:n
        if (SN(i).E<=0) && (SN(i).cond==1)  % if cluster heads energy depletes with reception
            SN(i).cond=0;
            SN(i).rop=rnd;
            dead_nodes=dead_nodes +1;
            operating_nodes= operating_nodes - 1;
            %             fprintf('IDL = %u \t , energy = %u \n', [i , SN(i).E].');
        end
        %              fprintf('IDNODE = %u \t , energy = %u \n', [i , SN(i).E].');
    end
    %     fprintf('====================================\n');
    
    if operating_nodes<n && temp_val==0
        temp_val=1;
        flag1stdead=rnd;
    end
    % Display Number of Cluster Heads of this round %
    %CLheads;
    
    
    transmissions=transmissions+1;
    if CLheads==0
        transmissions=transmissions-1;
    end
    
    
    % Next Round %
    rnd= rnd +1;
    fprintf('round = %u \n', [rnd].');
    tr(transmissions)=operating_nodes;
    op(rnd)=operating_nodes;
    
    if energy>0
        nrg(transmissions)=energy;
    end
    
    
    
    
    %     fprintf("OP = %u \n",[operating_nodes].');
    %     for i=1:n
    %         fprintf('ID = %u \t , energy = %u \n', [(i) , SN(i).E].');
    %     end
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
fprintf('OP = %u',[operating_nodes].');
% Plotting Simulation Results "Operating Nodes per Round" %
figure(4)
plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
title ({'LP'; 'Operating Nodes per Round';})
xlabel 'Rounds';
ylabel 'Operational Nodes';
hold on;

% Plotting Simulation Results  %
figure(5)
plot(1:transmissions,tr(1:transmissions),'-r','Linewidth',2);
title ({'LP'; 'Operational Nodes per Transmission';})
xlabel 'Transmissions';
ylabel 'Operational Nodes';
hold on;

% Plotting Simulation Results  %
figure(6)
plot(1:flag1stdead,nrg(1:flag1stdead),'-r','Linewidth',2);
title ({'LP'; 'Energy consumed per Transmission';})
xlabel 'Transmission';
ylabel 'Energy ( J )';
hold on;
