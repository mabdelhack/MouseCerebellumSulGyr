local_min = zeros(1,length(PCs));
local_max = zeros(1,length(PCs));
for i = 1 : length (PCs)
    local_min(i) = min(PCs{i}.Z);
    local_max(i) = max(PCs{i}.Z);
end
global_min = min(local_min);
global_max = max(local_max);

intervals = linspace(global_min,global_max,length(PCs));


for i = 1: length(PCs)
    quadrilaterals = DelaunayTri(PCs{i}.X,PCs{i}.Y,PCs{i}.Z);
    z_values = PCs{i}.Z(quadrilaterals.Triangulation);
    z_values_mean = mean(z_values,2);
    z_values_std = std(z_values,1,2);
    [n bins] = histc(z_values_mean,intervals);
    for j = 1 : length(intervals)-1
        slice = quadrilaterals.Triangulation(find ( bins == j & z_values_std<=1 ), : );
        union_shape = [];
        for k = 1 : size(slice,1)
            quadri.x = PCs{i}.X(slice(k,:));
            quadri.y = PCs{i}.Y(slice(k,:));
            triangles  = DelaunayTri(quadri.x,quadri.y);
            real_triangles(1).x =  quadri.x(triangles(1,:));
            real_triangles(1).y =  quadri.y(triangles(1,:));
            real_triangles(1).hole =  0;
            real_triangles(2).x =  quadri.x(triangles(2,:));
            real_triangles(2).y =  quadri.y(triangles(2,:));
            real_triangles(2).hole =  0;
            union_quad = PolygonClip(real_triangles(1),real_triangles(2),3);
            if k>1 
                union_shape = PolygonClip(union_shape,union_quad,3);
            else
                union_shape = union_quad;
            end
        end
        if isempty(union_shape)
            PCs_slices{i}{j} = [];
        else
            PCs_slices{i}{j} = PolygonClip(PCflattened{i},union_shape,1);
        end
        %%%%%%%%%%%%%%%%% testing code by plotting
%         if ~isempty(union_shape)
%             
%        
%             P1 = union_shape;
%             P2 = PCflattened{i};
%             figure
%             for type = 0:3
%                 subplot(2,2,type+1); box on
%                 switch type
%                     case 0; title('A-B')
%                     case 1; title('A.and.B (standard)')
%                     case 2; title('xor(A,B)')
%                     case 3; title('union(A,B)')
%                 end
%             P3=PolygonClip(P1,P2,type);
% 
%                 for io=1:3
%                     eval(['p=P' num2str(io) ';'])
%                     for np=1:length(p)
%                         obj=patch(p(np).x,p(np).y,io);
%                         if p(np).hole==1; set(obj,'facecolor','w'); end
%                     end
%                 end
% 
%             end
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of testing code
    end
end


%%%%%%%%%%%%%%%%%%%%%%% Code for computing the intersection
intersections_slice = zeros(length(PCs),length(PCs));
intersections_slice_normalized = zeros(length(PCs),length(PCs));
% intersection_out = cell(1,length(intervals));
for i = 1: length(PCs)-1
    for j = i+1:length(PCs)
        for k = 1: length(intervals)-1
            if (~isempty(PCs_slices{i}{k}) && ~isempty(PCs_slices{j}{k}))
                S(1).P = PCs_slices{i}{k};
                S(2).P = PCs_slices{j}{k};
%                 intersection_out{k} = Polygons_intersection(S,1,1e-6);
                intersection_out = Polygons_intersection(S,1,1e-6);
                if length(intersection_out) > 2
                    intersection_out_sec = intersection_out(3).area;
                else
                    intersection_out_sec = 0;
                end
            else
                intersection_out_sec = 0;
            end
            intersections_slice(i,j) = intersections_slice(i,j) + intersection_out_sec;    
        end
        normalizing_factor = (PCflatarea(i) + PCflatarea(j))/2;
        intersections_slice_normalized(i,j) = intersections_slice(i,j) / normalizing_factor; 
    end
end

Age = s1{1}(2:end-1);
Sulcus = strfind(s2{1},'sulcus');
if ~isempty(Sulcus)
    Region = 1;
else
    Region = 2;
end
PC1 = [];
PC2 = [];
Area = [];
NormalizedArea = [];
for i = 1:length(PCs)
    for j = i+1:length(PCs)
        PC1 = [PC1 i];
        PC2 = [PC2 j];
        Area = [Area; intersections_slice(i,j)];
        NormalizedArea = [NormalizedArea; intersections_slice_normalized(i,j)]
    end
end
Age = str2num(Age)*ones(length(PC1),1);
Position = Region*ones(length(PC1),1);
d = [PC1', PC2', Position, Age, Area, NormalizedArea];
dc = mat2cell(d,ones(1,size(d,1)),ones(1,size(d,2)));
dc = vertcat({'PC1', 'PC2', 'Position', 'Age', 'Intersection Area', 'Normalized Intersection Area'},dc);
cd       ../data_mtrs
cd       (s1{counter});
xlswrite([ s2{1}(9:end-3) 'xls'], dc, 'Intersection areas', 'A1');

PC = 1:length(PCs);
selfoverlapcount = selfoverlapcount';
selfoverlaparea  = selfoverlaparea';
PCflatarea = PCflatarea';
Age = Age(1)*ones(length(PC),1);
Position = Region*ones(length(PC),1);
d = [PC', Position, Age , PCflatarea, selfoverlapcount, selfoverlaparea  ];
dc = mat2cell(d,ones(1,size(d,1)),ones(1,size(d,2)));
dc = vertcat({'PC', 'Position', 'Age', 'Cell Area' 'Self-Intersection Count', 'Self-Intersection Area'},dc);
xlswrite([ s2{1}(9:end-3) 'self.xls'], dc, 'Intersection areas', 'A1');