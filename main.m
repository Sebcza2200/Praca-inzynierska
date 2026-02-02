clear
clc
close all

% UWAGA main.m z kodem źródłowym musi znajdować się w 
% tym samym folderze co pliki z danymi wejściowymi analizowanego samolotu 
% (i tylko analizowanego)

filename = 'dane.txt';
data = load(filename);

% UWAGA Proszę brać pod uwagę tylko takie wysokości w charakterystykach zespołu
% napędowego, przy których samolot uzyskuje przynajmniej jeden dodatni kąt
% toru lotu gama, w przeciwnym razie, wykres 3D nie zostanie wygenerowany, a
% jedynie wykresy 2D

fi = deg2rad(data(1));
a_inf=data(2);
S=data(3);
lambda_e=data(4);
a_0=deg2rad(data(5));
Cx_min=data(6);
n=data(7);
m_p=data(8);
m_0=data(9);
Cz_max=data(10);
q_e=data(11);

 ro_0=1.2255;
 g=9.81;
 ad_0=340.3;
 T_0=288.15;

disp('Możliwe do wyboru napędy samolotów');
disp('Samolot z napędem śmigłowym: 1');
disp('Samolot z napędem odrzutowym: 2');
wybor=input('Wybierz rodzaj napędu samolotu:  ');

dt=input('Wpisz wartość kroku czasowego dt [w sekundach]: ');
% Wartości 300 - 1000 są optymalne, 1000 - 3000 dla zgrubnych obliczeń

if wybor == 2
        
    files = dir('*.txt'); % Pobranie listy plików tekstowych w folderze
       
    % Posortowanie plików numerycznie      
    [~, sorted_idx] = sort(str2double(regexp({files.name}, '\d+', 'match', 'once')));
    files = files(sorted_idx);
    
    h_values = [];    
    Ma_all = [];    
    Pr_all = [];
          
    for k = 1:length(files)                
        filename = files(k).name; % Nazwa pliku    
               
        % Odczytanie danych z pliku               
        data = readlines(filename);    
              
        % Wysokość h (pierwsza linia pliku)             
        h = str2double(data(1));              
        h_values = [h_values; h];                 
                      
        measurements = data(2:end); % Odczytanie danych liczb Macha i ciągu (pozostałe linie)
                 
        Ma_temp = [];            
        Pr_temp = [];                    
               
        for i = 1:length(measurements)                     
                  
            % Rozdzielenie wartości na liczbę Macha i ciąg               
            line = strtrim(measurements(i)); % Usunięcie białych znaków        
                      
            if contains(line, ';')                          
                values = split(line, ';');                          
                if length(values) == 2                             
                    Ma_temp = [Ma_temp, str2double(values{1})]; % Liczba Macha                            
                    Pr_temp = [Pr_temp, str2double(values{2})]; % Ciąg                           
                else                             
                    warning('Nieprawidłowy format linii w pliku %s: %s', filename, line);                          
                end                       
            end              
        end    
             
        % Dopasowanie rozmiarów macierzy, uzupełnienie zerami              
        sizeV(k)=size(Ma_temp, 2);              
        max_measurements = max([size(Ma_all, 2), length(Ma_temp)]);              
        Ma_temp = [Ma_temp, zeros(1, max_measurements - length(Ma_temp))];             
        Pr_temp = [Pr_temp, zeros(1, max_measurements - length(Pr_temp))];    
               
        % Rozszerzenie wyników o nowe wiersze                 
        if isempty(Ma_all)                  
            Ma_all = Ma_temp;               
            Pr_all = Pr_temp;               
        else                 
            Ma_all = [Ma_all; zeros(1, size(Ma_all, 2))];                 
            Pr_all = [Pr_all; zeros(1, size(Pr_all, 2))];                 
            Ma_all(end, 1:length(Ma_temp)) = Ma_temp;                
            Pr_all(end, 1:length(Pr_temp)) = Pr_temp;              
        end           
    end

    V_all = zeros(length(h_values)-1, max(sizeV)); % Zainicjowanie macierzy V_all z odpowiednią liczbą wierszy i kolumn
    gama = zeros(length(h_values)-1, max(sizeV)); 
    w = zeros(length(h_values)-1, max(sizeV)); 

    for k = 1:length(files)-1
        h=h_values(k);
        tk=(m_p/n/q_e/(max(Pr_all(k,:))/1000));        
       
        for m = 1:tk/dt+1
            t = (m - 1) * dt;

            if h <= 11000            
                ro_h = ro_0 * (1 - h_values(k)/44300)^4.256;           
                T_h = T_0 - 6.5 * h_values(k)/1000;          
                ad_h = ad_0 * sqrt(T_h / T_0);       
            else           
                ro_h = 0.3637 * exp((h_values(k) - 11000) / 6337);          
                T_h = 216.65;          
                ad_h = 295.07;      
            end      
 
            V_minh(k, m) = sqrt(2 * g * (m_0 - t * n * q_e * max(Pr_all(k,:)) / 1000) / ro_h / S / Cz_max);
            i = sizeV(k); % Maksymalny rozmiar dla bieżącego k        
       
            while i > 0 && (Ma_all(k, i) * ad_h) > V_minh(k, m)
           
                if i <= size(Ma_all, 2) % Sprawdzenie zakresu               
                    V_all(k, i) = Ma_all(k, i) * ad_h; 

                    options = optimoptions('fsolve', 'Display', 'none');               
                    eqs = @(vars) [                   
                    0.5 * ro_h * S * a_inf * (vars(1) - a_0) * (Ma_all(k, i) * ad_h)^2 + n * Pr_all(k, i) * sin(vars(1) - fi) - g * (m_0 - t * n * q_e * (Pr_all(k, i) / 1000)) * cos(vars(2));                   
                    -0.5 * ro_h * S * ((Cx_min + 1 / lambda_e / 3.14 * a_inf^2 * (vars(1) - a_0)^2) / sqrt(1 - Ma_all(k, i)^2)) * (Ma_all(k, i) * ad_h)^2 + n * Pr_all(k, i) * cos(vars(1) - fi) - g * (m_0 - t * n * q_e * (Pr_all(k, i) / 1000)) * sin(vars(2))               
                    ];              
                    initial_guess = [0.01, 0.01];               
                    [R, ~, exitflag] = fsolve(eqs, initial_guess, options);               
                
                    if exitflag > 0                                            
                   
                        if isreal(R(1)) && isreal(R(2)) % Sprawdzenie, czy rozwiązanie jest zespolone                      
                            gama(k, i) = rad2deg(R(2));                       
                            w(k, i) = V_all(k, i) * sin(deg2rad(gama(k, i)));                   
                        else 
                            gama(k, i) = 0;                   
                        end               
                    else
                        gama(k, i) = 0;               
                    end           
                end
                
                % zapisanie wartości zmiennych w czasie
                V_all1(k,i+(m-1)*max(sizeV))=V_all(k,i);           
                gama_all(k,i+(m-1)*max(sizeV))=gama(k,i);
                w_all(k,i+(m-1)*max(sizeV))=w(k,i);            
                i = i - 1;        
        
            end         
        end 
    end

    wykresy2D(h_values, wybor, m_p, q_e, n, Pr_all, dt, sizeV, V_all1, gama_all, w_all);
    wykresyPowierzchniowe(h_values, wybor, m_p, q_e, n, Pr_all, dt, sizeV, V_all1, gama_all, w_all); 

elseif wybor == 1        
       
    files = dir('*.txt'); % Pobranie listy plików tekstowych w folderze
       
    % Posortowanie plików numerycznie      
    [~, sorted_idx] = sort(str2double(regexp({files.name}, '\d+', 'match', 'once')));      
    files = files(sorted_idx);
    
    h_values = [];      
    V_all = [];     
    Nr_all = [];
          
    for k = 1:length(files)
               
        filename = files(k).name; % Nazwa pliku  
        data = readlines(filename);  % Odczytanie danych z pliku   
            
        h = str2double(data(1));               
        h_values = [h_values; h];    
               
        % Odczytanie danych prędkości i mocy (pozostałe linie)              
        measurements = data(2:end);                
        V_temp = [];              
        Nr_temp = [];                    
               
        for i = 1:length(measurements)                     
                  
            % Rozdzielenie wartości na prędkości i moc                  
            line = strtrim(measurements(i)); % Usunięcie białych znaków        
                      
            if contains(line, ';')                          
                values = split(line, ';');                            
                if length(values) == 2                             
                    V_temp = [V_temp, str2double(values{1})]; % Prędkość                           
                    Nr_temp = [Nr_temp, str2double(values{2})]; % Moc                           
                else                             
                    warning('Nieprawidłowy format linii w pliku %s: %s', filename, line);                          
                end  
            end
        end
               
        % Dopasowanie rozmiarów macierzy, uzupełnienie zerami              
        sizeV(k)=size(V_temp, 2);              
        max_measurements = max([size(V_all, 2), length(V_temp)]);              
        V_temp = [V_temp, zeros(1, max_measurements - length(V_temp))];               
        Nr_temp = [Nr_temp, zeros(1, max_measurements - length(Nr_temp))];    
              
        % Rozszerzenie wyników o nowe wiersze  
        if isempty(V_all)                   
            V_all = V_temp;                
            Nr_all = Nr_temp;               
        else                   
            V_all = [V_all; zeros(1, size(V_all, 2))];                  
            Nr_all = [Nr_all; zeros(1, size(Nr_all, 2))];                  
            V_all(end, 1:length(V_temp)) = V_temp;                  
            Nr_all(end, 1:length(Nr_temp)) = Nr_temp;             
        end 
    end

    gama = zeros(length(h_values)-1, max(sizeV)); % Zainicjowanie macierzy
    w = zeros(length(h_values)-1, max(sizeV)); % Zainicjowanie macierzy

    for k = 1:length(files)-1
        h=h_values(k);
        tk=(m_p/n/(q_e/3600)/(max(Nr_all(k,:))/1000)); 
        
        for m = 1:1:tk/dt+1
            t = (m - 1) * dt;
        
            if h <= 11000           
                ro_h = ro_0 * (1 - h_values(k)/44300)^4.256;           
                T_h = T_0 - 6.5 * h_values(k)/1000;           
                ad_h = ad_0 * sqrt(T_h / T_0);       
            else          
                ro_h = 0.3637 * exp((h_values(k) - 11000) / 6337);           
                T_h = 216.65;           
                ad_h = 295.07;       
            end
       
            V_minh(k, m) = sqrt(2 * g * (m_0 - t * n * q_e / 3600 * max(Nr_all(k,:)) / 1000) / ro_h / S / Cz_max);       
            i = sizeV(k); % Maksymalny rozmiar dla bieżącego k        
       
            while i > 0 && V_all(k, i) > V_minh(k, m)
          
                if i <= size(V_all, 2) % Sprawdzenie zakresu
               
                    options = optimoptions('fsolve', 'Display', 'none');               
                    eqs = @(vars) [
                    0.5 * ro_h * S * a_inf * (vars(1) - a_0) * V_all(k, i)^2 + n * Nr_all(k, i) / V_all(k, i) * sin(vars(1) - fi) - g * (m_0 - t * n * (q_e / 3600) * (Nr_all(k, i) / 1000)) * cos(vars(2));                
                    -0.5 * ro_h * S * (Cx_min + 1 / lambda_e / 3.14 * a_inf^2 * (vars(1) - a_0)^2) * V_all(k, i)^2 + n * Nr_all(k, i) / V_all(k, i) * cos(vars(1) - fi) - g * (m_0 - t * n * (q_e / 3600) * Nr_all(k, i) / 1000) * sin(vars(2))
                    ];               
                    initial_guess = [0.01, 0.01];               
                    [R, ~, exitflag] = fsolve(eqs, initial_guess, options);
               
                    if exitflag > 0 
                        if isreal(R(1)) && isreal(R(2))   % Sprawdzenie, czy rozwiązanie jest zespolone  
                            gama(k, i) = rad2deg(R(2)); 
                            w(k, i) = V_all(k, i) * sin(deg2rad(gama(k, i)));  
                        else   
                            gama(k, i) = 0; 
                        end     
                    else     
                        gama(k, i) = 0; 
                    end   
                end
           
                V_all1(k,i+(m-1)*max(sizeV))=V_all(k,i);
                gama_all(k,i+(m-1)*max(sizeV))=gama(k,i);           
                w_all(k,i+(m-1)*max(sizeV))=w(k,i);           
                i = i - 1; 
            end 
        end
    end

    wykresy2D(h_values, wybor, m_p, q_e, n, Nr_all, dt, sizeV, V_all1, gama_all, w_all);
    wykresyPowierzchniowe(h_values, wybor, m_p, q_e, n, Nr_all, dt, sizeV, V_all1, gama_all, w_all);

else   
    disp("Nieprawidłowy wybór");
end

function wykresy2D(h_values, wybor, m_p, q_e, n, param1, dt, sizeV, V_all1, gama_all, w_all)

% Wykres funkcji gama(V) i w(V) dla różnych wysokości i czasów
for k = 1:length(h_values)-1
    all_negative = true; % Flaga do sprawdzania wartości gama_row
    first_plot = true; % Flaga do tworzenia figury tylko raz

    if wybor == 2
        tk = (m_p / n / q_e / (max(param1(k, :) / 1000)));
    else
        tk = (m_p / n / (q_e / 3600) / (max(param1(k, :) / 1000)));
    end

    colors = hsv(floor(tk / dt + 1)); % Generowanie unikalnych kolorów

    for m = 1:1:floor(tk / dt + 1)
        t = (m - 1) * dt; 
        mp_all(k, m) = m_p * (1 - t / tk);

        % Wyznaczenie zakresu indeksów
        start_idx = 1 + (m - 1) * max(sizeV);
        end_idx = min(start_idx + max(sizeV) - 1, size(V_all1, 2)); % Zabezpieczenie przed przekroczeniem rozmiaru
        if start_idx > size(V_all1, 2)
            break; % Przerwij pętlę, jeśli indeksy wykraczają poza rozmiar tablicy
        end

        % Pobranie danych
        V_row = V_all1(k, start_idx:end_idx);
        gama_row = gama_all(k, start_idx:end_idx);
        w_row = w_all(k, start_idx:end_idx);

        % Usunięcie zerowych wartości
        valid_idx = V_row > 0 & gama_row > 0;
        V_row = V_row(valid_idx);
        gama_row = gama_row(valid_idx);
        w_row = w_row(valid_idx); % Upewnij się, że w_row ma ten sam filtr

        % Sprawdzenie, czy są jakiekolwiek dodatnie wartości w gama_row
        if any(gama_row > 0)
            all_negative = false; % Przynajmniej jedna wartość jest dodatnia

            % Tworzenie figury i subplotów tylko raz
            if first_plot
                figure;
                subplot(2, 1, 1); % Pierwszy subplot na wykres gamma(V)
                hold on;
                first_plot = false;
            end

            % Rysowanie wykresu gamma(V) w pierwszym subplocie
            subplot(2, 1, 1);
            plot(V_row, gama_row, '-', 'Color', colors(m, :), 'DisplayName', sprintf('m_p = %.1f kg', mp_all(k, m)));

            % Rysowanie wykresu w(V) w drugim subplocie
            subplot(2, 1, 2);
            hold on;
            plot(V_row, w_row, '--', 'Color', colors(m, :), 'DisplayName', sprintf('m_p = %.1f kg', mp_all(k, m)));
        end           
    end
    
    if all_negative % Jeśli wszystkie wartości gama_row były ujemne, przerwij główną pętlę   
        break;
    end

    % Dodanie legendy, tytułu i opisów osi tylko, jeśli wykres został narysowany
    if ~first_plot
        % Ustawienia dla pierwszego subplotu (gamma(V))
        subplot(2, 1, 1);
        hold off;
        legend('show');
        xlabel('Prędkość V [m/s]');
        ylabel('Kąt drogi gamma [stopnie]');
        title(sprintf('Wykres gamma(V) dla różnych mas paliwa, na wysokości h = %.1f', h_values(k)));
        grid on;

        % Ustawienia dla drugiego subplotu (w(V))
        subplot(2, 1, 2);
        hold off;
        legend('show');
        xlabel('Prędkość V [m/s]');
        ylabel('Prędkość wznoszenia [m/s]');
        title(sprintf('Wykres w(V) dla różnych mas paliwa, na wysokości h = %.1f', h_values(k)));
        grid on;
    end
end
end

function wykresyPowierzchniowe(h_values, wybor, m_p, q_e, n, param1, dt, sizeV, V_all1, gama_all, w_all)
figure;
subplot(1, 2, 1); % Okno 1 - gamma
subplot(1, 2, 2); % Okno 2 - w

colors = hsv(length(h_values) - 1); % Kolory z palety tęczy

legend_entries = cell(1, length(h_values) - 1); % Komórka na teksty legendy

% Tworzymy uchwyty do powierzchni, żeby dodać legendę później
surface_handles_gamma = zeros(1, length(h_values) - 1);
surface_handles_w = zeros(1, length(h_values) - 1);

for k = 1:length(h_values)-1
    all_negative = true; % Flaga do sprawdzania wartości gama_row

    if wybor == 2
        tk = (m_p / n / q_e / (max(param1(k, :) / 1000)));
    else
        tk = (m_p / n / (q_e / 3600) / (max(param1(k, :) / 1000)));
    end

    % Przygotowanie danych do interpolacji
    V_all = [];
    gama_all_interp = [];
    mp_all_interp = [];
    w_all_interp = [];

    for m = 1:1:floor(tk / dt + 1)
        t = (m - 1) * dt; 
        mp_all(k, m) = m_p * (1 - t / tk);

        % Wyznaczenie zakresu indeksów
        start_idx = 1 + (m - 1) * max(sizeV);
        end_idx = min(start_idx + max(sizeV) - 1, size(V_all1, 2)); % Zabezpieczenie przed przekroczeniem rozmiaru
        if start_idx > size(V_all1, 2)
            break; % Przerwij pętlę, jeśli indeksy wykraczają poza rozmiar tablicy
        end

        % Pobranie danych
        V_row = V_all1(k, start_idx:end_idx);
        gama_row = gama_all(k, start_idx:end_idx);
        w_row = w_all(k, start_idx:end_idx);

        % Usunięcie zerowych wartości
        valid_idx = V_row > 0 & gama_row > 0;
        V_row = V_row(valid_idx);
        gama_row = gama_row(valid_idx);
        w_row = w_row(valid_idx); % Upewnij się, że w_row ma ten sam filtr

        % Sprawdzenie, czy są jakiekolwiek dodatnie wartości w gama_row
        if any(gama_row > 0)
            all_negative = false; % Przynajmniej jedna wartość jest dodatnia

            % Przechowywanie danych do interpolacji
            V_all = [V_all, V_row];
            gama_all_interp = [gama_all_interp, gama_row];
            mp_all_interp = [mp_all_interp, repmat(mp_all(k, m), size(V_row))]; % Replikacja masy dla każdego V_row
            w_all_interp = [w_all_interp, w_row]; % Replikacja wartości w dla każdego V_row
        end
    end

    if all_negative     % Jeśli wszystkie wartości gama_row były ujemne, przerwij główną pętlę
        break;
    end

    % Utworzenie siatki dla interpolacji w przedziale masy paliwa [0, m_p]
    [Xq, Yq] = meshgrid(linspace(min(V_all), max(V_all), 100), linspace(0, m_p, 100)); % Siatka dla interpolacji (masy paliwa od 0 do m_p)
    
    % Interpolacja gamma i w w zależności od masy paliwa i prędkości
    Zq_gamma = griddata(V_all, mp_all_interp, gama_all_interp, Xq, Yq, 'cubic'); % Interpolacja gamma
    Zq_w = griddata(V_all, mp_all_interp, w_all_interp, Xq, Yq, 'cubic'); % Interpolacja w

    % Wykres gamma(V,m) w pierwszym oknie (dla różnych wysokości)
    subplot(1, 2, 1);
    surface_handles_gamma(k) = surf(Xq, Yq, Zq_gamma, 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'FaceColor', colors(k, :));
    hold on;

    % Interpolacja punktów dla m_p na stałej wysokości
    idx_m_p = find(mp_all_interp == m_p); % Znalezienie punktów dla m_p
    V_m_p = V_all(idx_m_p); % Wartości prędkości dla m_p
    gama_m_p = gama_all_interp(idx_m_p); % Wartości gamma dla m_p
    plot3(V_m_p, repmat(m_p, size(V_m_p)), gama_m_p, '--', 'LineWidth', 2, 'Color','k'); % Łączenie punktów linią

    % Wykres w(V,m) w drugim oknie (dla różnych wysokości)
    subplot(1, 2, 2);
    surface_handles_w(k) = surf(Xq, Yq, Zq_w, 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'FaceColor', colors(k, :));
    hold on;

    % Interpolacja punktów dla m_p na stałej wysokości
    w_m_p = w_all_interp(idx_m_p); % Wartości w dla m_p
    plot3(V_m_p, repmat(m_p, size(V_m_p)), w_m_p, '--', 'LineWidth', 2, 'Color','k'); % Łączenie punktów linią

    % Dodanie wpisu do legendy (teraz używamy num2str)
    legend_entries{k} = [sprintf('h = %.1f m',h_values(k))]; % Upewnij się, że jest to tekst
end

% Dodanie etykiet i tytułów dla wykresu gamma
subplot(1, 2, 1);
xlabel('Prędkość V [m/s]');
ylabel('Masa paliwa m_p [kg]');
zlabel('Kąt drogi gamma [stopnie]');
title('Powierzchnia gamma(V,m) dla różnych wysokości');
legend(surface_handles_gamma, legend_entries, 'Location', 'Best'); % Dodanie legendy z wysokościami
zlim([0 inf]); % Zakres osi Z, zaczyna się od 0
grid on;

% Dodanie etykiet i tytułów dla wykresu w
subplot(1, 2, 2);
xlabel('Prędkość V [m/s]');
ylabel('Masa paliwa m_p [kg]');
zlabel('Prędkość wznoszenia [m/s]');
title('Powierzchnia w(V,m) dla różnych wysokości');
legend(surface_handles_w, legend_entries, 'Location', 'Best'); % Dodanie legendy z wysokościami
zlim([0 inf]); % Zakres osi Z, zaczyna się od 0
grid on;
end