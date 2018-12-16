%% this code provides us a quantitative measure for the amount of overlap neighboring Purkinje cells experience, 
% relative to precise anatomical region (sulcus vs.apex) and as a function of cellular growth (development).  
clear    all;
id       = [123];
s1       = {'P1501'};
area = 'sulcus';
s2       = { ...
      ['grouped_' s1{1} '_' area '.mtr']} % here we will add the rest of the files once they are traces

rotation = {[0 0 40]};  % rotated it so that the trees point up towards the molecular layer. 
scaling  = [80/40];  %shrinkage factor, slice thickness/ the shrunken thickness

counter  = 1; 
%% load the trees
cd       ./data_mtrs
cd       (s1{counter});
PCs      = load_tree (s2{counter});
PCs      = PCs {1}; % ungrouping the forest

load(['cut_point_' area '.mat'],'cut_point');

load(['P1_' area '.mat'],'P1');
load(['P2_' area '.mat'],'P2');
cd       ../../
%% delete somas to better see the arbors (ONLY USE FOR VISUALIZATION OTHERWISE THE OTHER SCRIPTS WILL NOT NOT NOT NOT NOT WORK !!!!!
% %   We don't need this turned on, as there are not too many points in the somata of these cells~~~~~~~
% for ward = 1 : length (PCs)
%     PCs{ward} = delete_tree (PCs{ward}, find (PCs{ward}.R ~= 1));  % does
%     not equal 1 will delete the 'tree' trace which in this case is
%     actually the Biological somata. 
% end 
for i = 1:length(PCs)
%     h1 = figure;
%     movegui(h1,'east');
%     xplore_tree(PCs{i});
%     h2 = figure;
%     movegui(h2,'west');
%     plot(PCs{i}.X,PCs{i}.Y);
    end_of_tree = length(PCs{i}.X);
    if ~isempty(cut_point(i))
        PCs{i} = delete_tree(PCs{i}, cut_point(i)+1:end_of_tree);
    end

%     close(h1);
%     close(h2);
end
%% translate forest, rotate it, scale it and plot it . !!!! Need to also resample it to have the same number of points!
clf;
T        = [PCs{1}.X(1) PCs{1}.Y(1) PCs{1}.Z(1)];
% rotation of whole forest around point T:
for ward       =  1 : length (PCs);
    
    % here need to resample the contours% ATTENTION 1
    
    PCs{ward}  = tran_tree  (PCs{ward}, -T); % translation T -> (0, 0, 0)
    PCs{ward}  = rot_tree   (PCs{ward},  rotation {counter});  % rotation
    % scales the trees in the z-plane:
    PCs{ward}  = scale_tree (PCs{ward}, [1 1 scaling(counter)]); 
    PCs{ward} = resample_tree(PCs{ward},4);
    PCs{ward}.D(:) = 1;
    hold on; 
    plot_tree(PCs{ward}, rand(1,3));   
end

%% ATTENTION 2: We have to fix the thickness of the contour/ tree to the 
% thickness. We can either do it in Neurolucida before we export the .xml
% file or we can automatize it here if the TREES toolbox has a function. 
for i = 1:length(PCs)
    PCs{i}.D(:) = 0;
    
end
% as you prefer and as you recommend.

%% ATTENTION 3: I think here we have to fix the intracellular crossing problem that you 
% and worked on already
% for i = 1:length(PCs)
%     h1 = figure;
%     movegui(h1,'east');
%     xplore_tree(PCs{i});
%     h2 = figure;
%     movegui(h2,'west');
%     plot(PCs{i}.X,PCs{i}.Y);
%  
% end
PCflatarea = zeros(1,length(PCs));
selfoverlaparea = zeros(1,length(PCs));
selfoverlapcount = zeros(1,length(PCs));
for i = 1:length(PCs)  
    
    U = [];
    for k = 1: length(P1{i})
        if size(P1{i}{k},1) > 1
            Pc1.x = [PCs{i}.X(P1{i}{k}(1,1):P1{i}{k}(1,2)) ; PCs{i}.X(P1{i}{k}(2,1):P1{i}{k}(2,2)); PCs{i}.X(P1{i}{k}(1,1))];
            Pc1.y = [PCs{i}.Y(P1{i}{k}(1,1):P1{i}{k}(1,2)) ; PCs{i}.Y(P1{i}{k}(2,1):P1{i}{k}(2,2)); PCs{i}.Y(P1{i}{k}(1,1))];
        else
            Pc1.x = PCs{i}.X(P1{i}{k}(1):P1{i}{k}(2));
            Pc1.y = PCs{i}.Y(P1{i}{k}(1):P1{i}{k}(2));
        end
        
        if size(P2{i}{k},1) > 1
            Pc2.x = [PCs{i}.X(P2{i}{k}(1,1):P2{i}{k}(1,2)) ; PCs{i}.X(P2{i}{k}(2,1):P2{i}{k}(2,2)); PCs{i}.X(P2{i}{k}(1,1))];
            Pc2.y = [PCs{i}.Y(P2{i}{k}(1,1):P2{i}{k}(1,2)) ; PCs{i}.Y(P2{i}{k}(2,1):P2{i}{k}(2,2)); PCs{i}.Y(P2{i}{k}(1,1))];
        else
            Pc2.x = PCs{i}.X(P2{i}{k}(1):P2{i}{k}(2));
            Pc2.y = PCs{i}.Y(P2{i}{k}(1):P2{i}{k}(2));
        end
        Intersection = PolygonClip(Pc1,Pc2,1);
        U  = [ U Intersection];
        selfoverlapcount(i) = selfoverlapcount(i) + length(Intersection);
        for op = 1: length(Intersection)
            selfoverlaparea(i) = selfoverlaparea(i) + polyarea(Intersection(op).x,Intersection(op).y);
        end
    end
    PCflattened{i}.x =[ PCs{i}.X ; PCs{i}.X(1)];
    PCflattened{i}.y =[ PCs{i}.Y ; PCs{i}.Y(1)];
    PCflattened{i}.hole = 0;
    if ~isempty(U)
        PCflattened{i} =PolygonClip(PCflattened{i},U,3);
    end
    
    for op = 1 : size(PCflattened{i},2)
        piece_area = (PCflattened{i}(op).hole * -2 +1 ) * polyarea(PCflattened{i}(op).x,PCflattened{i}(op).y);
        PCflatarea(i) = PCflatarea(i) + piece_area;
    end
    PCflatarea(i) = PCflatarea(i) + selfoverlaparea(i);
    figure;
    for np=1:length(PCflattened{i})
         obj=patch(PCflattened{i}(np).x,PCflattened{i}(np).y,'green');
         if PCflattened{i}(np).hole==1; set(obj,'facecolor','w'); end
    end
    for j = 1: length(U)
            patch(U(j).x,U(j).y,'red');
    end
    
end



