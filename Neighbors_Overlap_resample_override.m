%% this code provides us a quantitative measure for the amount of overlap neighboring Purkinje cells experience, 
% relative to precise anatomical region (sulcus vs.apex) and as a function of cellular growth (development).  
clear    all;
id       = [123];
s1       = {'P1502'};
area = 'apex';
s2       = { ...
      ['grouped_' s1{1} '_' area '_resampled.mtr']} % here we will add the rest of the files once they are traces

rotation = {[0 0 40]};  % rotated it so that the trees point up towards the molecular layer. 
scaling  = [80/40];  %shrinkage factor, slice thickness/ the shrunken thickness

counter  = 1; 
%% load the trees
cd       ../data_mtrs
cd       (s1{counter});
PCs      = load_tree (s2{counter});
% PCs      = PCs {1}; % ungrouping the forest

load(['cut_point_' area '.mat'],'cut_point');

load(['P1_' area '.mat'],'P1');
load(['P2_' area '.mat'],'P2');
cd       ../../Step2_Neighbor_Overlap
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

%     end_of_tree = length(PCs{i}.X);
%     if ~isempty(cut_point(i))
%         PCs{i} = delete_tree(PCs{i}, cut_point(i)+1:end_of_tree);
%     end

%     close(h1);
%     close(h2);
end
%% translate forest, rotate it, scale it and plot it . !!!! Need to also resample it to have the same number of points!
% clf;
% T        = [PCs{1}.X(1) PCs{1}.Y(1) PCs{1}.Z(1)];
% % rotation of whole forest around point T:
% for ward       =  1 : length (PCs);
%     
%     % here need to resample the contours% ATTENTION 1
%     
%     PCs{ward}  = tran_tree  (PCs{ward}, -T); % translation T -> (0, 0, 0)
%     PCs{ward}  = rot_tree   (PCs{ward},  rotation {counter});  % rotation
%     % scales the trees in the z-plane:
%     PCs{ward}  = scale_tree (PCs{ward}, [1 1 scaling(counter)]); 
%     PCs{ward} = resample_tree(PCs{ward},4);
%     hold on; 
%     plot_tree(PCs{ward}, rand(1,3));   
% end

%% ATTENTION 2: We have to fix the thickness of the contour/ tree to the 
% thickness. We can either do it in Neurolucida before we export the .xml
% file or we can automatize it here if the TREES toolbox has a function. 
% for i = 1:length(PCs)
%     PCs{i}.D(:) = 0;
% end
% % as you prefer and as you recommend.

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
        
%         %%%debugging
%         if i == 7
%             
%             figure;
%             hold on;
%             plot(Pc1.x,Pc1.y,'-r','linewidth',2);
%             plot(Pc2.x,Pc2.y,'-b','linewidth',2);
%         end
%         %%%end of debugging
        
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
%     figure;
%     for np=1:length(PCflattened{i})
%          obj=patch(PCflattened{i}(np).x,PCflattened{i}(np).y,'green');
%          if PCflattened{i}(np).hole==1; set(obj,'facecolor','w'); end
%     end
%     for j = 1: length(U)
%             patch(U(j).x,U(j).y,'red');
%     end
%     
end

%% ATTENTION 4: I think here is where we actually compute the overlap from a pairwise 
%type of method going throught the cells (I can suggest the method, but I
%have a hard time scripting, that's where your expertise comes in)

% overlap = cell(1,length(PCs));
% for i = 1: length(PCs)-1
%     for j = i+1:length(PCs)
%         if max(PCs{i}.Z) > min(PCs{j}.Z) | min(PCs{i}.Z) < max(PCs{j}.Z)
%             overlap{i} = [ overlap{i} j ];
%         end
%     end
% end
% 
% intersections = zeros(14,14);
% 
% for i = 1: length(PCs)-1
%     for j = i+1:length(PCs)
%         if max(PCs{i}.Z) > min(PCs{j}.Z) | min(PCs{i}.Z) < max(PCs{j}.Z)
%             clear S;
%             intersection_area = 0;
%             S(1).P = PCflattened{i};
%             
%             if (i == 1 && (j == 2 || j == 4))
%                 a = vertcat(PCflattened{i}.x, PCflattened{j}.x);
%                 maxx = max(a); minx = min(a); midx = (maxx + minx) / 2; midmidx = (maxx + midx) / 2;
%                 b = vertcat(PCflattened{i}.y, PCflattened{j}.y);
%                 maxy = max(b); miny = min(b); midy = (maxy + miny) / 2;
%                 S(2).P.x = [minx midx midx minx minx]; S(2).P.y = [miny miny midy midy miny]; S(2).P.hole = 0;
%                 inter_result = Polygons_intersection(S,0,1e-6);
%                 clear S;
%                 S(1).P = inter_result(3).P;
%                 S(2).P = PCflattened{j};
%                 intersection_out_sec = Polygons_intersection(S,0,1e-6);
%                 intersections(i,j) = intersections(i,j) + intersection_out_sec(3).area;
%                 clear S;
%                 S(1).P = inter_result(1).P;
%             end
%             S(2).P = PCflattened{j};
% %             S(4).P.x = [midx midmidx midmidx midx midx]; S(4).P.y = [miny miny midy midy miny]; S(4).P.hole = 0;
% %             S(5).P.x = [minx midx midx minx minx]; S(5).P.y = [midy midy maxy maxy midy]; S(5).P.hole = 0;
% %             S(6).P.x = [midx midmidx midmidx midx midx]; S(6).P.y = [midy midy maxy maxy midy]; S(6).P.hole = 0;
% %             S(7).P.x = [midmidx maxx maxx midmidx midmidx]; S(7).P.y = [miny miny midy midy miny]; S(7).P.hole = 0;
% %             S(8).P.x = [midmidx maxx maxx midmidx midmidx]; S(8).P.y = [midy midy maxy maxy midy]; S(8).P.hole = 0;
%             intersection_out =  Polygons_intersection(S,0,1e-6);
%             if length(intersection_out) > 2
%                 intersections(i,j) = intersections(i,j) + intersection_out(3).area;
%             end
%         end
%     end
% end
%% ATTENTION 5: I think here we will convert the output overlapping measure into a .txt 
%file arranged in a universal table so that other's 
% can understand clearly what we did, in the event that 
%anyone ever needs to see these results. (I can try to help here as well)

%% ATTENTION 6: here we save the results and table (I can give you some help here)





