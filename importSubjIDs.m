function subjlistBest10 = importSubjIDs(filename, dataLines)
%IMPORTFILE Import data from a text file
%  SUBJLISTBEST10 = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a string
%  array.
%
%  SUBJLISTBEST10 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  subjlistBest10 = importfile("/media/zhark2/GAOLAB_MRI/PipelineCmp/Toolbox_tmp/subjlist_Best10.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 23-Nov-2019 01:33:59

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = "BB1151_V2_26Feb2019";
opts.VariableTypes = "string";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "BB1151_V2_26Feb2019", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "BB1151_V2_26Feb2019", "EmptyFieldRule", "auto");

% Import the data
subjlistBest10 = readmatrix(filename, opts);

end