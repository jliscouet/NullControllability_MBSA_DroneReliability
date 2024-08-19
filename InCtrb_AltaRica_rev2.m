% =========================================================================
% Title: Integrate Controllability into AltaRica Models
% Authors: Zahra Motahari, Jonathan Liscouët
% Affiliation: Department of Mechanical, Industrial, and Aerospace Engineering
%              Concordia University, Montreal, Canada
% 
% Description:
% This MATLAB script integrates controllability into an existing 
% AltaRica model. Once the script is run, a GUI will prompt the user to:
% 1. Select a design from the design database.
% 2. Enter the maximum number of simultaneous failures to consider.
% 3. Select the path of the AltaRica model to use for integrating the 
%    controllability assessment.
% At the end of the process, the user will be asked to enter the name of 
% the updated AltaRica model file.
% 
% Inputs (via GUI):
% - selected_design (structure): The selected design structure.
% - n_f (double): The maximum number of simultaneous failures to consider.
% - selected_model_path (string): The file path to the selected AltaRica model.
% 
% Inputs (via command window at the end of the process):
% - outputModelName (string): The desired name for the output AltaRica 
%                             model, which will include the updated 
%                             controllability logics.
% 
% Outputs:
% - An AltaRica model file with the specified name, containing the 
%   integrated controllability logics.
% 
% Usage:
% 1. Ensure that the MATLAB environment is properly configured.
% 2. Run the script and follow the GUI prompts to select a design, specify 
%    the number of maximum simultaneous failures, and select the source 
%    AltaRica model file.
% 3. Enter the desired name for the output AltaRica model when prompted.
% 4. The script will generate the updated AltaRica model.
% 
% Example:
% [design, model_path, n_f] = selector();
% A = generate_controllable_cases(design, n_f);
% integrate_controllability(A, model_path);
% 
% Configuration Control:
% Version: 2.0
% Date: [2024-07-17]
% Change Log:
% - Version 1.0: Initial release of the script.
% - Version 2.0: Updated SELECTOR and generate_controllable_cases functions
% to include an option to select whether to include or exclude Yaw axis
% controllability.
% 
% License:
% This code is provided "as-is" under an open-source license. Users are 
% free to use, modify, and distribute this code with proper attribution to 
% the authors.
% =========================================================================

clear all; clc;

% Main script
[design, model_path,n_f,include_yaw] = selector();
A = generate_controllable_cases(design, n_f,include_yaw)
[num_rows, ~] = size(A);
disp(['Number of controllable cases: ', num2str(num_rows)]);
integrate_controllability(A,model_path)

%% ========================================================================
%% FUNCTION: selector()
%% ========================================================================
function [selected_design, selected_model_path, n_f, include_yaw] = selector()
    % SELECTOR Create a GUI for selecting a design and specifying the maximum simultaneous failures.
    %
    % This function creates a GUI that allows the user to select a design
    % by its name from the design database, stores the path of the selected AltArica model,
    % and specifies the number of maximum simultaneous failures.
    % Additionally, it allows the user to include or exclude Yaw axis controllability.

    % Load the designs from the DB_designs.m file
    designs = load_designs();
    
    % Extract the design names
    design_names = {designs.name};
    
    % Default number of maximum simultaneous failures
    default_n_f = 2;
    
    % Create the GUI figure
    fig = figure('Position', [300, 300, 500, 350], 'Name', 'Design Selector', 'NumberTitle', 'off', 'CloseRequestFcn', @close_gui);
    
    % Create a text label and popup menu for selecting a design
    uicontrol('Style', 'text', 'Position', [50, 260, 300, 30], 'HorizontalAlignment', 'left', 'String', '1. Select a Design:');
    popup = uicontrol('Style', 'popupmenu', 'Position', [350, 265, 130, 20], 'String', design_names);
    
    % Create a text label and input field for specifying the number of maximum simultaneous failures
    uicontrol('Style', 'text', 'Position', [50, 210, 300, 30], 'HorizontalAlignment', 'left', 'String', '2. Enter max simultaneous failures:');
    n_f_input = uicontrol('Style', 'edit', 'Position', [350, 215, 130, 20], 'String', num2str(default_n_f));
    
    % Create a text label and checkbox for Yaw axis controllability
    uicontrol('Style', 'text', 'Position', [50, 160, 300, 30], 'HorizontalAlignment', 'left', 'String', '3. Include Yaw Axis Controllability:');
    yaw_checkbox = uicontrol('Style', 'checkbox', 'Position', [350, 165, 130, 20], 'Value', 1);
    
    % Create a button to select a file
    uicontrol('Style', 'text', 'Position', [50, 110, 300, 30], 'HorizontalAlignment', 'left', 'String', '4. Select Model:');
    uicontrol('Style', 'pushbutton', 'Position', [350, 115, 130, 30], 'String', 'Select File', 'Callback', @select_file_callback);
    
    % Create a text label to display the selected model name
    model_name_label = uicontrol('Style', 'text', 'Position', [50, 80, 430, 20], 'HorizontalAlignment', 'left', 'String', 'Selected Model: ');
    
    % Create a button to store the selected model path and number of failures
    uicontrol('Style', 'pushbutton', 'Position', [200, 40, 100, 30], 'String', 'Done', 'Callback', @store_design_callback);
    
    % Variable to store the full path of the selected model
    selected_model_path = '';
    selected_design = [];
    n_f = default_n_f;
    include_yaw = true; % Default value for including Yaw axis controllability

    % Callback function for the Select File button
    function select_file_callback(~, ~)
        [file, path] = uigetfile('*.alt', 'Select AltArica File');
        if isequal(file, 0)
            disp('File selection canceled.');
        else
            selected_model_path = fullfile(path, file);
            [~, model_name, ext] = fileparts(selected_model_path);
            set(model_name_label, 'String', ['Selected Model: ', model_name, ext]);
        end
    end
    
    % Callback function for the Store Design button
    function store_design_callback(~, ~)
        % Get the selected design index
        selected_index = get(popup, 'Value');
        
        % Get the selected design
        selected_design = designs(selected_index);
        
        % Get the number of maximum simultaneous failures
        n_f = str2double(get(n_f_input, 'String'));
        
        % Validate the number of failures
        if isnan(n_f) || n_f <= 0
            errordlg('Please enter a valid number of maximum simultaneous failures.', 'Input Error');
            return;
        end
        
        % Get the Yaw axis controllability option
        include_yaw = get(yaw_checkbox, 'Value');
        
        % Check if a file has been selected
        if isempty(selected_model_path)
            errordlg('Please select a file.', 'File Selection Error');
            return;
        end
        
        % Display the selected design's details
        disp(['Selected Design: ', selected_design.name]);
        disp(['Selected Model: ', selected_model_path]);
        disp(['Number of Maximum Simultaneous Failures: ', num2str(n_f)]);
        disp(['Include Yaw Axis Controllability: ', num2str(include_yaw)]);
        
        % Resume the UI and close the GUI
        uiresume(fig);
        close(fig);
    end

    % Function to handle GUI close request
    function close_gui(~, ~)
        % Resume the UI in case it is closed before selection
        uiresume(fig);
        delete(fig);
    end

    % Wait for the user to make a selection
    uiwait(fig);
    
    % Output the selected design, model path, number of failures, and Yaw axis controllability option
    if isempty(selected_design) || isempty(selected_model_path)
        error('Selection incomplete.');
    end
end

function designs = load_designs()
    % LOAD_DESIGNS Load the designs from the DB_designs.m file
    run('DB_designs.m');
end

%% ========================================================================
%% FUNCTION: generate_controllable_cases()
%% ========================================================================
function A = generate_controllable_cases(design, n_f, include_yaw)
    % GENERATE_CONTROLLABLE_CASES Generate the matrix of controllable cases.
    %
    % Inputs:
    %   design - Structure containing B_f, u_min, u_max, and other parameters.
    %   n_f - Number of maximum simultaneous failures.
    %
    % Output:
    %   A - Matrix of controllable cases.

    % Extract parameters from the design structure
    B_f = design.B_f;
    u_min = design.u_min;
    u_max = design.u_max;
    G = -[-design.m * 9.789, 0, 0, 0]'; % Disturbance vector, in [T L M N]'

    % Remove Yaw Axis Consideration if Requested
    if include_yaw == 0
        B_f(end, :) = [];  % Remove B_f yaw row (last row)
        G(end) = [];       % Remove G yaw component (last component)
    end

    % Generate Failure Matrix
    [~, num_columns] = size(B_f);
    W1 = generate_combinations(num_columns, n_f);
    
    % Initialize variables
    [num_rows, ~] = size(W1);
    A = [];

    for counter = 1:num_rows
        W = diag(W1(counter, :));
        
        % Initialize Control effectiveness matrix with failure case
        B_f_failed = B_f * W;

        % ACAI computation
        delta = u_max - u_min; % Linear effector contribution magnitude
        fc = (u_max + u_min) / 2; % Define center of envelope
        Fc = B_f_failed * fc; % Center point of envelope
        [n_Bf, m_Bf] = size(B_f_failed); % Dimension of B_f_failed
        S1 = nchoosek(1:m_Bf, n_Bf - 1); % Combinations of (n-1) effectors out of M

        % Compute ACAI for each hyperplane segment
        d = zeros(1, size(S1, 1));
        for i = 1:size(S1, 1)
            choose = S1(i, :); % Combination of (n-1) effectors
            B1 = B_f_failed(:, choose); % B1 matrix (hypersegment row-space)
            xi = null(B1'); % xi vector (orth. to hypersegment col-space)
            xi = xi(:, 1) / norm(xi(:, 1)); % Normalized xi vector
            B2f = B_f_failed;
            B2f(:, choose) = [];
            delta2 = delta;
            delta2(choose) = [];
            l = abs(xi' * B2f) * abs(delta2); % Sum of projected absolute contribution
            g = abs(xi' * (Fc - G)); % Absolute projected distance from G to Fc
            d(i) = l / 2 - g; % Difference between G-to-Fc and hyperplane segment-to-Fc
        end

        % Determine if the failure case is controllable
        rho = min(d);
        if rho > 1e-6
            A = [A; diag(W)'];
        end
    end
end

%% ========================================================================
%% FUNCTION: generate_combinations()
%% ========================================================================
function W1 = generate_combinations(n, k)
    % GENERATE_COMBINATIONS Generate combinations of n elements with at most k zeros.
    % This function generates a matrix W1 where each row is a unique
    % combination of n elements (each being 0 or 1) with the constraint
    % that the combination has at most k elements being 0.
    %
    % Inputs:
    %   n - Number of elements in each combination (number of columns).
    %   k - Maximum number of elements being 0 in each combination.
    %
    % Outputs:
    %   W1 - A matrix containing all valid combinations.

    % Initialize an empty cell array to store valid combinations
    valid_combinations = {};
    
    % Generate combinations using nested loops
    for num_zeros = 0:k
        % Generate all combinations with exactly num_zeros zeros
        comb = nchoosek(1:n, num_zeros);
        for i = 1:size(comb, 1)
            combination = ones(1, n);
            combination(comb(i, :)) = 0;
            valid_combinations{end + 1} = combination;
        end
    end
    
    % Convert the cell array to a matrix
    W1 = cell2mat(valid_combinations');
end

%% ========================================================================
%% FUNCTION: integrate_controllability()
%% ========================================================================

function integrate_controllability(controllable_cases_matrix, model_path)
    % INTEGRATE_CONTROLLABILITY Integrate controllable cases into an AltArica model.
    %
    % Inputs:
    %   controllable_cases_matrix - Matrix of controllable cases.
    %   model_path - Path to the source AltArica model file.

M = controllable_cases_matrix;
sourceModelPath = model_path;

% 0. Prompt user to enter the name of the output model file
[~, sourceFileName, ~] = fileparts(sourceModelPath);
defaultOutputName = [sourceFileName, 'ctrb.alt'];
prompt = 'Enter the name for the output AltaRica model file: ';
outputModelName = input([prompt '(' defaultOutputName '): '], 's');

% Use the default name if the user provides no input
if isempty(outputModelName)
    outputModelName = defaultOutputName;
end

% 1. Access number of controllable cases and rotors
[L, R] = size(M);                                                           % L is the number of controllable cases (-), R is the rotor number (-) 

% 2. Access Altarica code
% Altarica Model Path
fid = fopen(sourceModelPath,'rt') ;                                         % Load source file
X = fread(fid) ;                                                            % X is the full code
fclose(fid) ;                                                               % Close code file
X = char(X.') ;                                                             % Convert X from [text] to [charact]

% 3. Locate start of controllability block code
A = strfind ( X ,'block Controllability');                                  % Find block title
A1 = strfind ( X ,'observer'); % Find the end of code
B = X( A + 26 : A1 );                                                        % Start editing 26 characters after block title
A2 = strfind (B,'end'); % Find the end of controllabilty block
A2 = A2(1); 
B = B(1:A2); % Extract text of the controllability block
A3 = strfind (B,'Integer Out'); % Find the Integer Out in the text of controllability block
B2 = B(A3:A2); % Extract "Integer Out_XXX (reset = 0); ↵ assertion ↵ Out_XXX := 0;"
A4 = strfind(B2,';');
A4 = A4(1);
B2 = B2 (1:A4); % Extract the "Integer Out_XXX (reset = 0);" from the text of controllability block
B = erase(B,B2); % Remove the "Integer Out_XXX (reset = 0);" from the text of controllability block
A5 = strfind(B, 'assertion');
B = B(1:A5-1);
A6 = strfind (B2,'Out'); 
A7 = strfind (B2,'(');
B3 = B2(A6:A7-1); % Extract "Out_XXX"
A8 = strfind (B, 'Integer') + 8; % Find the number of inputs of controllability block
A9 = strfind (B, '(') -1; % A8 and A9 find the interval of Input number, for exampler: Integer _1_192 (reset=0) --> find _1_192


% 4. Access Controllability block input names
CC = {};                                                                    % CC is set of input names separated by ","
for i = 1 : 1 : R
    if i < R
    C = B (A8(i): A9(i));
    CC{2*(i-1)+1} = char(C);
    CC{2*(i-1)+2} = ',';
    else
    C = B (A8(i): A9(i));
    CC{2*(i-1)+1} = char(C);
    end
end

% 5. Assign a logic integer to each controllable case
for i = 1:L
    if i ==1
        str = sprintf('  Integer Logic_%d (reset = 0); ',i);                % str is the text declaring a logic integer
    else
        str = [str sprintf(' Integer Logic_%d (reset = 0); ',i)];
    end
end

% 6. ???
k = 1;                                                                      % K is counter for align the position of variable is S cell array
for i = 1:L
    for j = 1:R
        if M(i,j) == 1
            S{i,k} = CC(2*j-1);                                             % S is the cell array of the sets of the input flow variables involved in each controllable failure case
            % Odd positions belong to the input name
            if j < R
                AA = sum (M(i,j+1:R)) ;
                if AA > .5
                
            S{i,k+1} = CC(2*j);
                else
                end
            % Even positions belong to the commas after input name
            k = k+2;
            else
            end
        else
        end
    end
    k = 1;
end

% 7. Create the controllable case assertions
Z = textscan(str,'%s');                                                     % Z is a cell array with the data in str (logic integer declarations)
for i = 1:L  
    SS{i} = [Z{1}{5*(i-1)+2} ':= max (' S{i,:} ');' ];                      % SS contains the controllable case assertions. For example, Logic_1 := max(_1_192, ..., _1_198); 
end

% 8. Integrate the controllable case assertions with AND operator
f = round(sqrt(L));                                                         % Limit the number of inputs to the test assertions, because of combinatorial explosion of tests to run for compilation

if L > 3.5

for i = 1:f                                                                 
F = '  := min (Test%d, ';
F1 = [B3,F];
      if i == 1
         str1 = sprintf(F1,i); 
      else
          if i < f
         str1 = [str1 sprintf(' Test%d, ',i)];
          else
             str1 = [str1 sprintf(' Test%d ',i)];          
          end

      end

end

 

for i = 1:f % Make the Logic of each scenario that it defines with Integers in Altarica

    if i ==1
        str2 = sprintf('Integer  Test%d , ',i);
    elseif i<f
        str2 = [str2 sprintf(' Test%d , ',i)];
    else
        str2 = [str2 sprintf(' Test%d (reset = 0); ',i)];
    end

end

Z1 = textscan (str2,'%s');                                                  % Convert cell to text
DD = SS;


for counter =1:f

   
    for i = f*(counter-1)+1 :1:f*counter                                    % Fill the Test by Logics e.g Test1 := min (Logic_1, ..., Logic_8);
    if i < L && i < f*counter 
        DDs{counter,i,i+1} = [DD{1,i}(1),DD{1,i}(4)];
         SSS{counter} = [Z1{1}(2*counter) ':= min (' DDs{counter,:} ');' ];
    elseif i == L || i == f*counter && i < L+0.5
 DDs{counter,i,i+1} = [DD{1,i}(1),[]];
        SSS{counter} = [Z1{1}(2*counter) ':= min (' DDs{counter,:} ');' ];
    else 
        break
    end

    end

end

elseif L>1.5
    SSS = [];
    for i = 1:L                                                                 
F = '  := min (Logic_%d, ';
F1 = [B3,F];
      if i == 1
         str1 = sprintf(F1,i); 
      else
          if i < L
         str1 = [str1 sprintf(' Logic_%d, ',i)];
          else
             str1 = [str1 sprintf(' Logic_%d ',i)];          
          end

      end

end


    str2 = [];
    f = 0;
else
    F = '  :=  (Logic_%d ';
F1 = [B3,F];
    str1 = sprintf(F1,1); 
    str2 = [];
    f = 0;
    SSS = [];
end
SS = [SS SSS];

% 9. ? Bring together each part of controllability block assertion
S22 = {};
str1 = [str1 sprintf('); ')];                                               
S2 = [str newline str2 newline 'assertion  ' newline SS{1:L+f} newline str1];            % Integration
S22 = join (S2);

% 10. Replace string S1 with string S2
F2 = ':= 0;';
F3 = [B3,F2];
V1 = ['assertion' newline ,'            ',F3];                             
Y = strrep(X, V1, S22) ;

% 11. Edit destination file
fid2 = fopen(outputModelName,'wt');                                         % Access destination file
fwrite(fid2,Y{1});                                                          % Update destination file
fclose(fid2);                                                               % Close destination file

% 12. Inform the user that the process is completed
[~, sourceFileName, sourceFileExt] = fileparts(sourceModelPath);
sourceFileNameWithExt = [sourceFileName, sourceFileExt];
disp(['The integration of Null Controllability into ', sourceFileNameWithExt, ' is complete.']);
disp(['The updated model has been saved as: ', outputModelName]);
end
