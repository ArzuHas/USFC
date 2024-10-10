clear;
clc;

% Constants
num_aal = 90; % replace with appropriate number of brain regions specific to the atlas image that is aimed to be used 
listpath = fullfile('path_to_lists_directory');  % Replace with appropriate path for the subject list
datapath = fullfile('path_to_fc_matrices');  % Replace with appropriate path for FC matrices
outputpath = fullfile('path_to_output_directory');  % Replace with appropriate path
% SC = structural connectivity
% FC = functional connectivity

% Load necessary data as .mat file
load('AAL_ind_116to90.mat'); % update index (ind): AAL atlas compressed to 90 investigated regions such removing cerebellum and vermis from the atlas

% Import subject IDs with both FC and SC matrices
subjects_both = importSubjIDs(fullfile(listpath, 'HCP_1200_list_both_fmask.txt')); %replace with appropriate file name for subject list
num_s = length(subjects_both);

% Load or compute center of mass and distance matrices
outfile1 = 'CMs_AAL_adult_90AAL.mat'; %replace with the atlas specific indices
if ~exist(outfile1, 'file')
    AAL_data = niftiread('AAL.nii'); % replace with appropriate atlas image that is aimed to be used 
    ROIs = unique(AAL_data);
    ROIs(ROIs == 0) = []; % remove background
    CMs = zeros(length(ROIs), 3);
    for idx = 1:length(ROIs)
        [x, y, z] = ind2sub(size(AAL_data), find(AAL_data == ROIs(idx)));
        CMs(idx, :) = [mean(x), mean(y), mean(z)];
    end
    CM_90 = CMs(ind, :);
    save(outfile1, 'CMs', 'CM_90', 'ind');
end

outfile2 = 'Distance_CMs_AAL_adult_90AAL.mat';
if ~exist(outfile2, 'file')
    load(outfile1);
    Distance_90 = zeros(num_aal, num_aal); %replace with appropriate number of brain regions for a specific atlas 
    for rr = 1:num_aal
        for rrr = rr:num_aal
            Distance_90(rr, rrr) = norm(CM_90(rr, :) - CM_90(rrr, :));
            Distance_90(rrr, rr) = Distance_90(rr, rrr);
        end
    end
    save(outfile2, 'Distance_90');
end

% Main computation loop
for ss = 1:num_s
    fprintf('Calculating for subject %d\n', ss);
    subj = subjects_both{ss};

    % Define file paths
    finfile = fullfile(outputpath, 'FC_Matrices_FDR_corrected_90AAL', [subj '_Functional_Matrix_FDR_corrected.mat']);%replace with appropriate file name for FC
    sinfile = fullfile(outputpath, 'Structural_Matrices_individual_nothr_90AAL', [subj '_Structural_Matrix.mat']); %%replace with appropriate file name for SC
    outputdir_USFC = fullfile('path_to_output_USFC_directory');  % Replace with appropriate path

    if ~exist(outputdir_USFC, 'dir')
        mkdir(outputdir_USFC);
    end

    % Check if output files already exist
    outfile_all = fullfile(outputdir_USFC, [subj '_Cost_Route_all_Matrix.mat']);
    if ~exist(outfile_all, 'file')
        load(finfile);
        load(sinfile);
     % Make the economical assumption to find the most efficient segments up to 4 steps
        Cost_M = Distance_90 ./ Structural_M;
        Min_cost_M = inf(num_aal, num_aal);
        Route_M = cell(num_aal, num_aal);

        for rr = 1:(num_aal - 1)
            for rrr = (rr + 1):num_aal
                if Functional_M(rr, rrr) ~= 0
                    % Find minimum cost routes (up to 4 steps)
                    min_routes = cell(4, 1);
                    min_costs = zeros(4, 1);
                    min_costs(1) = Cost_M(rr, rrr);
                    min_routes{1} = [rr, rrr];

                    % 2 steps
                    [min_cost, I1] = min(Cost_M(rr, :) + Cost_M(:, rrr)');
                    min_costs(2) = min_cost;
                    min_routes{2} = [rr, I1, rrr];

                    % 3 steps
                    costs = arrayfun(@(x, y) Cost_M(rr, x) + Cost_M(x, y) + Cost_M(y, rrr), 1:num_aal, 1:num_aal);
                    [min_cost, I] = min(costs(:));
                    [I1, I2] = ind2sub(size(costs), I);
                    min_costs(3) = min_cost;
                    min_routes{3} = [rr, I1, I2, rrr];

                    % 4 steps
                    costs = arrayfun(@(x, y, z) Cost_M(rr, x) + Cost_M(x, y) + Cost_M(y, z) + Cost_M(z, rrr), 1:num_aal, 1:num_aal, 1:num_aal);
                    [min_cost, I] = min(costs(:));
                    [I1, I2, I3] = ind2sub(size(costs), I);
                    min_costs(4) = min_cost;
                    min_routes{4} = [rr, I1, I2, I3, rrr];

                    % Select the minimum route
                    [All_min_cost, All_I] = min(min_costs);
                    Min_cost_M(rr, rrr) = All_min_cost;
                    if ~isinf(All_min_cost)
                        Route_M{rr, rrr} = min_routes{All_I};
                    end
                end
            end
        end

        % Symmetrize matrices
        It = logical(tril(ones(num_aal, num_aal), -1));
        Min_cost_M(It) = Min_cost_M';
        Route_M(It) = Route_M';

        % Save the results
        save(outfile_all, 'Cost_M', 'Min_cost_M', 'Route_M');
    end

    % Route-specific matrices
    Route_M1_eff = cell(num_aal, num_aal);
    Route_M2_eff = cell(num_aal, num_aal);
    Route_M3_eff = cell(num_aal, num_aal);
    Route_M4_eff = cell(num_aal, num_aal);

    for rr = 1:(num_aal-1)
        for rrr = rr+1:num_aal
            route_len = length(Route_M{rr, rrr});
            switch route_len
                case 2
                    Route_M1_eff{rr, rrr} = Route_M{rr, rrr};
                case 3
                    Route_M2_eff{rr, rrr} = Route_M{rr, rrr};
                case 4
                    Route_M3_eff{rr, rrr} = Route_M{rr, rrr};
                case 5
                    Route_M4_eff{rr, rrr} = Route_M{rr, rrr};
            end
        end
    end

    % Symmetrize the route matrices
    It = logical(tril(ones(num_aal,num_aal),-1));
    Route_M1_eff(It) = Route_M1_eff';
    Route_M2_eff(It) = Route_M2_eff';
    Route_M3_eff(It) = Route_M3_eff';
    Route_M4_eff(It) = Route_M4_eff';

    % Save route-specific matrices
    save(fullfile(outputdir_USFC, [subj '_Route1_Matrix.mat']), 'Route_M1_eff');
    save(fullfile(outputdir_USFC, [subj '_Route2_Matrix.mat']), 'Route_M2_eff');
    save(fullfile(outputdir_USFC, [subj '_Route3_Matrix.mat']), 'Route_M3_eff');
    save(fullfile(outputdir_USFC, [subj '_Route4_Matrix.mat']), 'Route_M4_eff');

    % Further processing and saving FC and SC matrices
    outfile_FC_M1 = fullfile(outputdir_USFC, [subj '_FC_R1_Matrix.mat']);
    outfile_FC_M2 = fullfile(outputdir_USFC, [subj '_FC_R2_Matrix.mat']);
    outfile_FC_M3 = fullfile(outputdir_USFC, [subj '_FC_R3_Matrix.mat']);
    outfile_FC_M4 = fullfile(outputdir_USFC, [subj '_FC_R4_Matrix.mat']);
    outfile_USFC_all = fullfile(outputdir_USFC, [subj '_USFC_Matrix.mat']);
    outfile_USFC_all_abs = fullfile(outputdir_USFC, [subj '_USFC_Matrix_abs.mat']);
    outfile_Route_count = fullfile(outputdir_USFC, [subj '_USFC_Route_count_Matrix.mat']);

    if ~exist(outfile_FC_M1, 'file')
        load(finfile);
        load(sinfile);

        SC_M1 = zeros(num_aal, num_aal);
        SC_M2 = zeros(num_aal, num_aal);
        SC_M3 = zeros(num_aal, num_aal);
        SC_M4 = zeros(num_aal, num_aal);  
        FC_M1 = zeros(num_aal, num_aal);
        FC_M2 = zeros(num_aal, num_aal);
        FC_M3 = zeros(num_aal, num_aal);
        FC_M4 = zeros(num_aal, num_aal);
       
        % Construct the SC and FC for each step 
        for rr = 1:num_aal
            for rrr = 1:num_aal
                Route = Route_M{rr, rrr};
                if ~isempty(Route)
                    n_steps = length(Route);
                    switch n_steps
                        case 2
                            SC_M1(rr, rrr) = Structural_M(Route(1), Route(2));
                            FC_M1(rr, rrr) = Functional_M(Route(1), Route(2));
                        case 3
                            SC_M2(rr, rrr) = (Structural_M(Route(1), Route(2)) + Structural_M(Route(2), Route(3))) / 2;
                            FC_M2(rr, rrr) = Functional_M(Route(1), Route(3));
                        case 4
                            SC_M3(rr, rrr) = (Structural_M(Route(1), Route(2)) + Structural_M(Route(2), Route(3)) + Structural_M(Route(3), Route(4))) / 3;
                            FC_M3(rr, rrr) = Functional_M(Route(1), Route(4));
                        case 5
                            SC_M4(rr, rrr) = (Structural_M(Route(1), Route(2)) + Structural_M(Route(2), Route(3)) + Structural_M(Route(3), Route(4)) + Structural_M(Route(4), Route(5))) / 4;
                            FC_M4(rr, rrr) = Functional_M(Route(1), Route(5));
                    end
                end
            end
        end

        % Symmetrize matrices
        SC_M1(It) = SC_M1';
        SC_M2(It) = SC_M2';
        SC_M3(It) = SC_M3';
        SC_M4(It) = SC_M4';
        FC_M1(It) = FC_M1';
        FC_M2(It) = FC_M2';
        FC_M3(It) = FC_M3';
        FC_M4(It) = FC_M4';

        % Save the results
        save(outfile_FC_M1, 'FC_M1');
        save(outfile_FC_M2, 'FC_M2');
        save(outfile_FC_M3, 'FC_M3');
        save(outfile_FC_M4, 'FC_M4');
    end

    % Calculate and save USFC matrices
    if ~exist(outfile_USFC_all, 'file')
        USFC_M = zeros(num_aal, num_aal);
        USFC_M_abs = zeros(num_aal, num_aal);
        RouteCounts_M = zeros(num_aal, num_aal);

        for rr = 1:num_aal
            for rrr = 1:num_aal
                Route = Route_M{rr, rrr};
                if ~isempty(Route)
                    n_steps = length(Route) - 1;
                    for nn = 1:n_steps
                        step_i = Route(nn);
                        step_j = Route(nn+1);
                        if step_i < step_j
                            USFC_M(step_i, step_j) = USFC_M(step_i, step_j) + Functional_M(rr, rrr);
                            USFC_M_abs(step_i, step_j) = USFC_M_abs(step_i, step_j) + abs(Functional_M(rr, rrr));
                            RouteCounts_M(step_i, step_j) = RouteCounts_M(step_i, step_j) + 1;
                        elseif step_j < step_i
                            USFC_M(step_j, step_i) = USFC_M(step_j, step_i) + Functional_M(rrr, rr);
                            USFC_M_abs(step_j, step_i) = USFC_M_abs(step_j, step_i) + abs(Functional_M(rrr, rr));
                            RouteCounts_M(step_j, step_i) = RouteCounts_M(step_j, step_i) + 1;
                        end
                    end
                end
            end
        end

        % Symmetrize matrices
        USFC_M(It) = USFC_M';
        USFC_M_abs(It) = USFC_M_abs';
        RouteCounts_M(It) = RouteCounts_M';

        % Save the results
        save(outfile_USFC_all, 'USFC_M');
        save(outfile_USFC_all_abs, 'USFC_M_abs');
        save(outfile_Route_count, 'RouteCounts_M');
    end

    %  Save the SC and FC matrices as edge files for further analysis and plotting (replace with appropriate paths)

    baseFolder_SC_M = fullfile('path_to_SC_M_directory'); 
    if ~exist(baseFolder_SC_M, 'dir')
        mkdir(baseFolder_SC_M);
    end
    saveAsEdge(SC_M1, fullfile(baseFolder_SC_M, [subj, '_SC_M1.edge']));
    saveAsEdge(SC_M2, fullfile(baseFolder_SC_M, [subj, '_SC_M2.edge']));
    saveAsEdge(SC_M3, fullfile(baseFolder_SC_M, [subj, '_SC_M3.edge']));
    saveAsEdge(SC_M4, fullfile(baseFolder_SC_M, [subj, '_SC_M4.edge']));

    baseFolder_FC_M = fullfile('path_to_FC_M_directory');
    if ~exist(baseFolder_FC_M, 'dir')
        mkdir(baseFolder_FC_M);
    end
    saveAsEdge(FC_M1, fullfile(baseFolder_FC_M, [subj, '_FC_M1.edge']));
    saveAsEdge(FC_M2, fullfile(baseFolder_FC_M, [subj, '_FC_M2.edge']));
    saveAsEdge(FC_M3, fullfile(baseFolder_FC_M, [subj, '_FC_M3.edge']));
    saveAsEdge(FC_M4, fullfile(baseFolder_FC_M, [subj, '_FC_M4.edge']));

    baseFolder_USFC_M = fullfile('path_to_USFC_M_directory'); 
    if ~exist(baseFolder_USFC_M, 'dir')
        mkdir(baseFolder_USFC_M);
    end
    saveAsEdge(USFC_M, fullfile(baseFolder_USFC_M, [subj, '_USFC.edge']));
    saveAsEdge(USFC_M_abs, fullfile(baseFolder_USFC_M, [subj, '_USFC_abs.edge']));

    baseFolder_Route = fullfile('path_to_Route_directory');
    if ~exist(baseFolder_Route, 'dir')
        mkdir(baseFolder_Route);
    end
    saveAsEdge(RouteCounts_M, fullfile(baseFolder_Route, [subj, '_Route_counted.edge']));
end
