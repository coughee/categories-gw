clear all;
myFolder = 'C:\cygwin\home\Jonathan\Projects\gwcatsandgods\';
filePattern = fullfile(myFolder, '*.txt');
txtFiles = dir(filePattern);

for k = 1:length(txtFiles)
    baseFileName = txtFiles(k).name;
    fullFileName = fullfile(myFolder,baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    godtemp = dlmread(fullFileName);
    
end