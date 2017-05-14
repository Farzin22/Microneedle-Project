%bargraph production code
%author Farzin Akbar

colormap(gray)
%the_path =strcat(data_path,'\','*.csv');
v=strfind(bar_path,'\')
the_path=bar_path(1:end-(length(bar_path)-v(length(v)))-1)
cd(the_path);
%dd = dir('*.csv');                           
dd = dir(bar_path);                           

%%%%%%%replacing comma with dots for the decimal point because German computer!
fileNames = {dd.name}; 
copyfile(dd.name,'temp_bar.csv')
for i=1:length(fileNames)

	file    = memmapfile( 'temp_bar.csv', 'writable', true );
	comma   = uint8(',');
	point   = uint8('.');
	file.Data( transpose( file.Data==comma) ) = point;

end

%%%%%%%%%%%%%%%%%%%%%%%%%saving the replaced values in another temporary%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = cell(numel(fileNames),2);
data(:,1) = regexprep(fileNames, '.csv','');

for ii = 1:numel(fileNames)    
   data{ii,2} = dlmread('temp_bar.csv',';',1,0);
   [n,s,r] = xlsread('temp_bar.csv');
end

f = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(f, 'Position'); % grab original on screen size
set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(f,'PaperPositionMode','auto') %set paper pos for printing
diameters =data{1,2}(:,1);
diameters_surrogate =[];
for i=1:length(diameters)
    diameters_surrogate(i)= i
end

%data{1,2}(:,1);

subplot(2,2,1) %penetration force bars
Penetration_Force_Average_mN = data{1,2}(:,2);
Penetration_Force_std_mN = data{1,2}(:,3);
bar(diameters_surrogate,Penetration_Force_Average_mN,'facecolor','r'); %drawing the bars at the bar graph
set(gca,'XTickLabel',diameters)
hold on 
errorbar(diameters_surrogate,Penetration_Force_Average_mN,Penetration_Force_std_mN,'.','Linewidth',2,'color','k') %drawing the error bars at the bar graph
xlabel('Tip Diameter (µm)');
ylabel('Penetration Force (mN)');
grid on
str1 = ['$\bar{F}_{penetration}$ =', sprintf('%5.2f mN, ', Penetration_Force_Average_mN), '$\sigma$=', sprintf('%5.2f mN, ', Penetration_Force_std_mN)];

subplot(2,2,2) %force drop bars
set(gca,'fontsize',20)
Force_Drop_Average_mN = data{1,2}(:,6);
Force_Drop_std_mN = data{1,2}(:,7);
bar(diameters_surrogate,Force_Drop_Average_mN,'facecolor','c'); %drawing the bars at the bar graph
set(gca,'XTickLabel',diameters)
hold on 
errorbar(diameters_surrogate,Force_Drop_Average_mN,Force_Drop_std_mN,'.','Linewidth',2,'color','k') %drawing the error bars at the bar graph
xlabel('Tip Diameter (µm)');
ylabel('Force Drop (mN)');
grid on
str2 = ['$\bar{F}_{drop}$ =', sprintf('%5.2f mN, ', Force_Drop_Average_mN), '$\sigma$=', sprintf('%5.2f mN, ', Force_Drop_std_mN)];

subplot(2,2,3) %displacement at penetration bars
set(gca,'fontsize',20)
Displacement_at_Penetration_Average = data{1,2}(:,10);
Displacement_at_Penetration_std = data{1,2}(:,11);
bar(diameters_surrogate,Displacement_at_Penetration_Average,'facecolor','m'); %drawing the bars at the bar graph
set(gca,'XTickLabel',diameters)
hold on 
errorbar(diameters_surrogate,Displacement_at_Penetration_Average,Displacement_at_Penetration_std,'.','Linewidth',2,'color','k') %drawing the error bars at the bar graph
xlabel('Tip Diameter (µm)');
ylabel('Displacement at Penetration(µm)');
grid on
str3 = ['$\bar{D}_{penetration}$ =', sprintf('%5.2f, ', Displacement_at_Penetration_Average),'$\mu$ ', '$ m   \sigma $ = ', sprintf('%5.2f, ', Displacement_at_Penetration_std)];

subplot(2,2,4) % time at penetration bars
set(gca,'fontsize',20)
Time_at_Penetration_Average = data{1,2}(:,14);
Time_at_Penetration_std = data{1,2}(:,15);
bar(diameters_surrogate,Time_at_Penetration_Average,'facecolor','g'); %drawing the bars at the bar graph
set(gca,'XTickLabel',diameters)
hold on
errorbar(diameters_surrogate,Time_at_Penetration_Average,Time_at_Penetration_std,'.','Linewidth',2,'color','k') %drawing the error bars at the bar graph
xlabel('Tip Diameter (µm)');
ylabel('Time at Penetration(s)');
grid on
str4 = ['$\bar{T}_{penetration}$ =', sprintf('%5.2f s, ', Time_at_Penetration_Average), '$ \sigma $ = ', sprintf('%5.2f s, ', Time_at_Penetration_std)];
%title(str4,'interpreter','latex')


subplot(2,2,1)
set(gca,'fontsize',20)
subplot(2,2,2)
set(gca,'fontsize',20)
subplot(2,2,3)
set(gca,'fontsize',20)
subplot(2,2,4)
set(gca,'fontsize',20)

file_name =strcat(bar_path(1:end-4),'\','bar_graphs');
new_folder=char(fileNames);
mkdir(new_folder(1:end-4))
cd(bar_path(1:end-4))
saveas(f, 'bar_graphs') % save figure
print(f,'bar_graphs','-dpng', '-r300')
print(f,'bar_graphs','-dtiff')
set(f,'Position', origSize) %set back to original dimensions
to_print_bargraph=strcat(bar_path(1:end-4),'\','bar_graphs.png');

fclose('all');

