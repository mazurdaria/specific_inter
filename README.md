
# Численное моделирование двойного электрического слоя на границе металл-электролит

## Обзор
Код предназначен для решения уравнений самосогласованного поля для расчета профилей концентраций ионов и электростатического потенциала в растворах электролитов и ионных жидкостях вблизи плоского бесконечного электрода с учетом структурных взаимодействий в рамках теории типа Кана-Хиллиарда и исключенного объема ионов в рамках моделей несимметричного решеточного газа или двухкомпонентной смеси твердых сфер в рамках приближения Перкуса-Йевика. При решении уравнений используется встроенная функция bvp4c MATLAB для решения краевых задач. Код разработан в среде MATLAB 2022a и не требует дополнительных установок

## Особенности
Решает уравнения химического равновесия ионов и уравнение Пуассона для потенциала электростатического поля в растворах электролитов и ионных жидкостях вблизи плоского бесконечного заряженного электрода. Уравнения самосогласованного поля получены с помощью вариации функционала термодинамического потенциала. Теория учитывает структурные взаимодействия ионов в рамках подхода Кана-Хиллиарда через билинейную форму по градиентам концентраций ионов в термодинамическом потенциале. Позволяет получать результаты для двух моделей: двухкомпонентной смеси твердых сфер и несимметричного решеточного газа 
Позволяет учитывать специфическую адсорбцию ионов на плоский бесконечный электрод через дополнительный потенциал внешнего поля.

## Начало работы

### Предварительные требования
- MATLAB 2022a

### Установка
Дополнительная установка не требуется. Просто склонируйте или загрузите этот репозиторий на свой локальный компьютер.

### Использование
Чтобы запустить решатель, откройте MATLAB скрипт в MATLAB 2022a и выполните его. Скрипт автоматически рассчитает профили концентрации и потенциала на основе предопределенных параметров и граничных условий.

```matlab
% Для запуска скрипта используйте командное окно MATLAB
run('path_to_script.m')
```

### Вывод
Скрипт выводит:
- Профили концентраций ионов вблизи электрода. 
- Профиль электростатического потенциала. 
- Графические диаграммы для визуализации этих профилей.

## Лицензия
Этот проект распространяется под [лицензией MIT](LICENSE.md).

## Литература
Разработанный код основан на следующих публикациях: Для обеих моделей решается система (30): https://www.researchgate.net/publication/372355956_Theory_of_electrolyte_solutions_in_a_slit_charged_pore_Effects_of_structural_interactions_and_specific_adsorption_of_ions Однако для случая модели смеси твердых сфер, химические потенциалы компонент имеют вид согласно уравнениям (45) и (46) в статье: https://www.researchgate.net/publication/360570219_Modified_Poisson-Boltzmann_equations_and_macroscopic_forces_in_inhomogeneous_ionic_fluids)

---
