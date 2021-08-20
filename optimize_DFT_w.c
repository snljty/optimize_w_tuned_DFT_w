/*******************************************************************************
 * optimize w, or literally omega, in long-range correction functional of DFT. *
 *******************************************************************************/

/* the application of Brent's method is part of github.com:rqtl/qtl2scan */

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# ifdef _WIN32
# include <io.h>
# else
# include <unistd.h>
# endif
# include <math.h>
# include <time.h>

# include "brent_fmin.h"

int glob_argc = 1;
long glob_replace_iop_pos = 0l;
char glob_gau_exe[BUFSIZ + 1] = "";

# define Close_file(flp) fclose(flp); flp = NULL

void Print_exit_success();
void Print_exit_failure();
void Pause_program(char const *prompt);
void Get_J_and_J_squared(double *J_ptr, double *J_squared_ptr);
double Calc_J_squared_from_w(double w, void *args);

int main(int argc, char const *argv[])
{
    unsigned int iarg = 0u;

    double w_low = 0.05, w_high = 0.6;
    double w_guess = (w_low + w_high) / 2;
    double w_tolerance = 1E-4;

    double w_when_J_squared_min = 0.0;

    unsigned int multi_n = 0u, multi_np1 = 0u, multi_nm1 = 0u; /* p for +, n for -, 0 for not set */
    int charge_n = 0, charge_np1 = -1, charge_nm1 = 1;

    unsigned int const max_iter = 100u;

    char const temp_name[] = "template.gjf";
    FILE *temp_ifl = NULL;
    FILE *n_ifl = NULL, *np1_ifl = NULL, *nm1_ifl = NULL;

    char env_gauss_exedir_copy[BUFSIZ + 1] = "";
    char *env_gauss_exedir = NULL;
    char const *env_gauss_exedir_ptr = getenv("GAUSS_EXEDIR");

    char buf[BUFSIZ + 1] = "";
    char *line_end = NULL;

    int info = 0;

    # ifdef _WIN32
    char const path_splitter[] = ";";
    # else
    char const path_splitter[] = ":";
    # endif

    time_t time_start = 0, time_stop = 0;

    glob_argc = argc;

    /* check command arguments, if "-h" or "--help" appears, print help message and exit. */
    for (iarg = 1; iarg != argc; ++ iarg)
    {
        if (! strcmp(argv[iarg], "-h") || ! strcmp(argv[iarg], "--help") || ! strcmp(argv[iarg], "/?"))
        {
            printf("Usage: %s [ARGUMENTS]\n", argv[0]);
            printf("ARGUMENTS: \n");
            printf("    [ -h | --help ]                         Print this help message and exit.\n");
            printf("    [ --low w_LOW ]                         The lower limit of w.\n");
            printf("    [ --high w_HIGH ]                       The higher limit of w.\n");
            printf("    [ --guess w_GUESS ]                     The initial guess of w.\n");
            printf("    [ --multi-np1 MULTIPLICITY_N+1 ]        The multiplicity of N+1 state.\n");
            printf("    [ --multi-nm1 MULTIPLICITY_N-1 ]        The multiplicity of N-1 state.\n");
            printf("    [ --tolerance TOLERANCE ]               The tolerance of convergence of w.\n");
            printf("\n");
            printf("\"N\" stands for the reference state, \"N+1\" stands for \"N\" plus an extra electron, \n");
            printf("and \"N-1\" stands for \"N\" minus an electron.\n");
            printf("The default is: W_LOW = %4.2lf, W_HIGH = %4.2lf\n", w_low, w_high);
            printf("W_GUESS = (W_LOW + W_HIGH) / 2, MULTIPLICITY_N-1 and MULTIPLICITY_N+1 = (both) \n");
            printf("\"multiplicity of reference state\" + 1, ");
            printf("where \"multiplicity of reference state\" is read from \"template.gjf\", and\n");
            printf("TOLERANCE = %6.1lg.\n", w_tolerance);
            printf("\n");
            printf("You need to prepare a template file called \"template.gjf\" in the current working directory, \n");
            printf("which is the entire input single point energy task file, except the IOps for tuning w.\n");
            printf("\n");
            Print_exit_success();
        }
    }

    /* check template file */
    temp_ifl = fopen(temp_name, "rt");
    if (! temp_ifl)
    {
        fprintf(stderr, "Error! Cannot find \"%s\".\n", temp_name);
        Print_exit_failure();
    }
    Close_file(temp_ifl);

    /* parse command arguments */
    iarg = 0;
    for (;;)
    {
        ++ iarg;
        if (iarg == argc)
            break;
        if (! strcmp(argv[iarg], "--low"))
        {
            ++ iarg;
            if (iarg == argc)
            {
                fprintf(stderr, "Error! Missing argument after \"%s\".\n", argv[iarg - 1]);
                Print_exit_failure();
            }
            if (sscanf(argv[iarg], "%lg", & w_low) != 1)
            {
                fprintf(stderr, "Error! Cannot recognize value \"%s\" after \"%s\".\n", argv[iarg], argv[iarg - 1]);
                Print_exit_failure();
            }
            if (w_low <= 0 || (unsigned int)(w_low * 1E4) < 1)
            {
                fprintf(stderr, "Error! Minimal acceptable w is 0.0001, but got %6.1lg.\n", w_low);
                Print_exit_failure();
            }
            w_guess = (w_low + w_high) / 2;
            continue;
        }
        if (! strcmp(argv[iarg], "--high"))
        {
            ++ iarg;
            if (iarg == argc)
            {
                fprintf(stderr, "Error! Missing argument after \"%s\".\n", argv[iarg - 1]);
                Print_exit_failure();
            }
            if (sscanf(argv[iarg], "%lg", & w_high) != 1)
            {
                fprintf(stderr, "Error! Cannot recognize value \"%s\" after \"%s\".\n", argv[iarg], argv[iarg - 1]);
                Print_exit_failure();
            }
            if (w_high <= 0.0)
            {
                fprintf(stderr, "Error! Higher limit of w must be positive, but got %6.1lg.\n", w_high);
                Print_exit_failure();
            }
            if ((unsigned int)(w_high * 1E4) >= 1E5)
            {
                fprintf(stderr, "Error! Maximum acceptable w is 10, but got %6.1lg.\n", w_high);
                Print_exit_failure();
            }
            w_guess = (w_low + w_high) / 2;
            continue;
        }
        if (! strcmp(argv[iarg], "--guess"))
        {
            ++ iarg;
            if (iarg == argc)
            {
                fprintf(stderr, "Error! Missing argument after \"%s\".\n", argv[iarg - 1]);
                Print_exit_failure();
            }
            if (sscanf(argv[iarg], "%lg", & w_guess) != 1)
            {
                fprintf(stderr, "Error! Cannot recognize value \"%s\" after \"%s\".\n", argv[iarg], argv[iarg - 1]);
                Print_exit_failure();
            }
            continue;
        }
        if (! strcmp(argv[iarg], "--multi-np1"))
        {
            ++ iarg;
            if (iarg == argc)
            {
                fprintf(stderr, "Error! Missing argument after \"%s\".\n", argv[iarg - 1]);
                Print_exit_failure();
            }
            if (sscanf(argv[iarg], "%u", & multi_np1) != 1)
            {
                fprintf(stderr, "Error! Cannot recognize value \"%s\" after \"%s\".\n", argv[iarg], argv[iarg - 1]);
                Print_exit_failure();
            }
            if (! multi_np1)
            {
                fprintf(stderr, "Error! Multiplicity cannot be zero, but it is zero for N+1 state.\n");
                Print_exit_failure();
            }
            continue;
        }
        if (! strcmp(argv[iarg], "--multi-nm1"))
        {
            ++ iarg;
            if (iarg == argc)
            {
                fprintf(stderr, "Error! Missing argument after \"%s\".\n", argv[iarg - 1]);
                Print_exit_failure();
            }
            if (sscanf(argv[iarg], "%u", & multi_nm1) != 1)
            {
                fprintf(stderr, "Error! Cannot recognize value \"%s\" after \"%s\".\n", argv[iarg], argv[iarg - 1]);
                Print_exit_failure();
            }
            if (! multi_nm1)
            {
                fprintf(stderr, "Error! Multiplicity cannot be zero, but it is zero for N-1 state.\n");
                Print_exit_failure();
            }
            continue;
        }
        if (! strcmp(argv[iarg], "--tolerance"))
        {
            ++ iarg;
            if (iarg == argc)
            {
                fprintf(stderr, "Error! Missing argument after \"%s\".\n", argv[iarg - 1]);
                Print_exit_failure();
            }
            if (sscanf(argv[iarg], "%lg", & w_tolerance) != 1)
            {
                fprintf(stderr, "Error! Cannot recognize value \"%s\" after \"%s\".\n", argv[iarg], argv[iarg - 1]);
                Print_exit_failure();
            }
            if (w_tolerance < 1E-4)
            {
                fprintf(stderr, "Error! Minimum acceptable tolerance of w is 0.0001, but got %6.1lg.\n", w_tolerance);
                Print_exit_failure();
            }
            continue;
        }
        fprintf(stderr, "Error! Cannot recognize argument \"%s\".\n", argv[iarg]);
        Print_exit_failure();
    }
    if (w_high <= w_low)
    {
        fprintf(stderr, "Error! higher limit of w (%6.1lg) must be greater than lower limit of w (%6.1lg).\n", w_high, w_low);
        Print_exit_failure();
    }
    if (w_guess < w_low || w_guess > w_high)
    {
        fprintf(stderr, "Error! Initial guess of w (%6.1lg) must be in interval [lower limit (%6.1lg), higher limit (%6.1lg)].\n", \
            w_guess, w_low, w_high);
        Print_exit_failure();
    }


    if (! env_gauss_exedir_ptr || ! strcmp(env_gauss_exedir_ptr, ""))
    {
        fprintf(stderr, "Error! Environment variable \"GAUSS_EXEDIR\" is not set properly!\n");
        Print_exit_failure();
    }
    strcpy(env_gauss_exedir_copy, env_gauss_exedir_ptr);
    env_gauss_exedir = strtok(env_gauss_exedir_copy, path_splitter);
    while (env_gauss_exedir)
    {
        if (strchr(env_gauss_exedir, ' '))
        {
            # ifdef _WIN32
            sprintf(glob_gau_exe, "\"%s%s\"", env_gauss_exedir, "\\g16.exe"); /* assume it is Gaussian 16 */
            # else
            sprintf(glob_gau_exe, "\"%s%s\"", env_gauss_exedir, "/g16");
            # endif
        }
        else
        {
            # ifdef _WIN32
            sprintf(glob_gau_exe, "%s%s", env_gauss_exedir, "\\g16.exe");
            # else
            sprintf(glob_gau_exe, "%s%s", env_gauss_exedir, "/g16");
            # endif
        }
        if (access(glob_gau_exe, X_OK))
        {
            if (strchr(env_gauss_exedir, ' '))
            {
                # ifdef _WIN32
                sprintf(glob_gau_exe, "\"%s%s\"", env_gauss_exedir, "\\g09.exe"); /* assume it is Gaussian 16 */
                # else
                sprintf(glob_gau_exe, "\"%s%s\"", env_gauss_exedir, "/g09");
                # endif
            }
            else
            {
                # ifdef _WIN32
                sprintf(glob_gau_exe, "%s%s", env_gauss_exedir, "\\g09.exe");
                # else
                sprintf(glob_gau_exe, "%s%s", env_gauss_exedir, "/g09");
                # endif
            }
            if (! access(glob_gau_exe, X_OK))
                break;
        }
        else
            break;
        env_gauss_exedir = strtok(NULL, path_splitter);
    }
    if (access(glob_gau_exe, X_OK))
    {
        fprintf(stderr, "Error! Cannot find either g16 or g09 as executable.\n");
        Print_exit_failure();
    }

    /* start preparing input files */
    temp_ifl = fopen(temp_name, "rt");
    n_ifl = fopen("N.gjf", "wt");
    np1_ifl = fopen("Np1.gjf", "wt");
    nm1_ifl = fopen("Nm1.gjf", "wt");
    /* link 0 ,route section and a blank line followed */
    while (fgets(buf, BUFSIZ, temp_ifl))
    {
        if ((line_end = strstr(buf, "\r\n")))
            strcpy(line_end, "\n");
        if (* buf == '#')
        {
            * strchr(buf, '\n') = '\0';
            fprintf(n_ifl, "%s", buf);
            fprintf(np1_ifl, "%s", buf);
            fprintf(nm1_ifl, "%s", buf);
            glob_replace_iop_pos = ftell(n_ifl);
            fprintf(n_ifl, "%s", " IOp(3/107=0000000000,3/108=0000000000)\n");
            fprintf(np1_ifl, "%s", " IOp(3/107=0000000000,3/108=0000000000)\n");
            fprintf(nm1_ifl, "%s", " IOp(3/107=0000000000,3/108=0000000000)\n");
            continue;
        }
        fprintf(n_ifl, "%s", buf);
        fprintf(np1_ifl, "%s", buf);
        fprintf(nm1_ifl, "%s", buf);
        if (* buf == '\n')
            break;
    }
    /* title and a blank line followed */
    while (fgets(buf, BUFSIZ, temp_ifl))
    {
        if ((line_end = strstr(buf, "\r\n")))
            strcpy(line_end, "\n");
        fprintf(n_ifl, "%s", buf);
        fprintf(np1_ifl, "%s", buf);
        fprintf(nm1_ifl, "%s", buf);
        if (* buf == '\n')
            break;
    }
    /* charge and multiplicity */
    fgets(buf, BUFSIZ, temp_ifl);
    if ((line_end = strstr(buf, "\r\n")))
        strcpy(line_end, "\n");
    if (sscanf(buf, "%d %u", & charge_n, & multi_n) != 2)
    {
        fprintf(stderr, "Error! Cannot read the charge and multiplicity from template file.\n");
        Close_file(temp_ifl);
        Close_file(n_ifl);
        Close_file(np1_ifl);
        Close_file(nm1_ifl);
        remove("N.gjf");
        remove("Np1.gjf");
        remove("Nm1.gjf");
        Print_exit_failure();
    }
    if (! multi_n)
    {
        fprintf(stderr, "Error! Multiplicity cannot be zero, but it is zero in template file for N state.\n");
        Close_file(temp_ifl);
        Close_file(n_ifl);
        Close_file(np1_ifl);
        Close_file(nm1_ifl);
        remove("N.gjf");
        remove("Np1.gjf");
        remove("Nm1.gjf");
        Print_exit_failure();
    }
    charge_np1 = charge_n - 1;
    charge_nm1 = charge_n + 1;
    if (! multi_np1)
        multi_np1 = multi_n + 1;
    else
    {
        if (! ((multi_np1 - multi_n) & 1))
        {
            fprintf(stderr, "Error! Multiplicity of N+1 state and N state must have different parity.\n");
            Close_file(temp_ifl);
            Close_file(n_ifl);
            Close_file(np1_ifl);
            Close_file(nm1_ifl);
            remove("N.gjf");
            remove("Np1.gjf");
            remove("Nm1.gjf");
            Print_exit_failure();
        }
    }
    if (! multi_nm1)
        multi_nm1 = multi_n + 1;
    else
    {
        if (! ((multi_nm1 - multi_n) & 1))
        {
            fprintf(stderr, "Error! Multiplicity of N-1 state and N state must have different parity.\n");
            Close_file(temp_ifl);
            Close_file(n_ifl);
            Close_file(np1_ifl);
            Close_file(nm1_ifl);
            remove("N.gjf");
            remove("Np1.gjf");
            remove("Nm1.gjf");
            Print_exit_failure();
        }
    }
    fprintf(n_ifl, "%d %u\n", charge_n, multi_n);
    fprintf(np1_ifl, "%d %u\n", charge_np1, multi_np1);
    fprintf(nm1_ifl, "%d %u\n", charge_nm1, multi_nm1);
    /* atom coordinates and others */
    while (fgets(buf, BUFSIZ, temp_ifl))
    {
        fprintf(n_ifl, "%s", buf);
        fprintf(np1_ifl, "%s", buf);
        fprintf(nm1_ifl, "%s", buf);
    }
    /* close these files */
    Close_file(temp_ifl);
    Close_file(n_ifl);
    Close_file(np1_ifl);
    Close_file(nm1_ifl);

    /* show title */
    printf("Optimize w (literally omega) in long-range correction functional of DFT.\n");
    printf("Parameters: w_low = %6.4lf, w_high = %6.4lf, w_guess = %6.4lf, w_tolerance = %6.4lf\n", \
        w_low, w_high, w_guess, w_tolerance);
    printf("            Charges for N, N+1 and N-1 states: %d %d %d\n", charge_n, charge_np1, charge_nm1);
    printf("            Multiplicities for N, N+1 and N-1 states: %u %u %u\n", multi_n, multi_np1, multi_nm1);
    printf("\n");
    time_start = time(NULL);

    /* Brent's method for minimize J^2 with variable w. */
    w_when_J_squared_min = Brent_fmin(w_low, w_high, w_guess, Calc_J_squared_from_w, NULL, w_tolerance, max_iter, & info);
    if (info > 0)
    {
        fprintf(stderr, "Error! Arguments of Brent's method are illegal!\n");
        Print_exit_failure();
    }
    else if (info < 0)
    {
        fprintf(stderr, "Error! w did not converge, the last value is: %6.4lf\n", w_when_J_squared_min);
        time_stop = time(NULL);
        printf("Total time elapsed: %d s.\n", (int)difftime(time_stop, time_start));
        exit(EXIT_FAILURE);
    }
    printf("Minimum value of J^2 encountered when w = %6.4lf.\n", w_when_J_squared_min);
    printf("You can use \"IOp(3/107=%05u00000,3/108=%05u00000)\" in your further Gaussian input files.\n", \
        (unsigned int)(w_when_J_squared_min * 1E4), (unsigned int)(w_when_J_squared_min * 1E4));
    printf("\n");

    /* end of task */
    time_stop = time(NULL);
    printf("Total time elapsed: %d s.\n", (int)difftime(time_stop, time_start));
    printf("\n");
    remove("N.gjf");
    remove("N.out");
    remove("Np1.gjf");
    remove("Np1.out");
    remove("Nm1.gjf");
    remove("Nm1.out");

    /* pause program on Windows is no command arguments are provided. */
    # ifdef _WIN32
    if (argc == 1)
        Pause_program("Press <Enter> to exit ...");
    # endif
    Print_exit_success();

    return 0;
}

void Print_exit_success()
{
    fprintf(stdout, "Exiting normally.\n");
    exit(EXIT_SUCCESS);
}

void Print_exit_failure()
{
    fprintf(stderr, "Exiting abnormally.\n");
    # ifdef _WIN32
    if (glob_argc == 1)
        Pause_program("Press <Enter> to exit ...");
    # endif
    exit(EXIT_FAILURE);

    return;
}

void Pause_program(char const *prompt)
{
    char cont = '\0';

    if (prompt)
        printf("%s\n", prompt);
    while ((cont = getchar()) != '\n' && cont != EOF)
        ;

    return;
}

void Get_J_and_J_squared(double *J_ptr, double *J_squared_ptr)
{
    double J_n = 0.0, J_np1 = 0.0;
    double E_n = 0.0, E_np1 = 0.0, E_nm1 = 0.0;
    double e_HOMO_n = 0.0, e_HOMO_np1 = 0.0;
    FILE *n_ofl = NULL, *np1_ofl = NULL, *nm1_ofl = NULL;
    long last_line_pos = 0l, this_line_pos = 0l;
    char last_value_str[BUFSIZ + 1] = "";
    char buf[BUFSIZ + 1] = "";
    char *tok = NULL;

    /* open files */
    n_ofl = fopen("N.out", "rt");
    if (! n_ofl)
    {
        fprintf(stderr, "Error! Cannot open \"N.out\" for reading! Check your Gaussian settings.");
        Print_exit_failure();
    }
    np1_ofl = fopen("Np1.out", "rt");
    if (! np1_ofl)
    {
        fprintf(stderr, "Error! Cannot open \"Np1.out\" for reading! Check your Gaussian settings.");
        Close_file(n_ofl);
        Print_exit_failure();
    }
    nm1_ofl = fopen("Nm1.out", "rt");
    if (! nm1_ofl)
    {
        fprintf(stderr, "Error! Cannot open \"Nm1.out\" for reading! Check your Gaussian settings.");
        Close_file(n_ofl);
        Close_file(np1_ofl);
        Print_exit_failure();
    }

    /* electron energy */
    /* electron energy of N */
    while (fgets(buf, BUFSIZ, n_ofl))
    {
        if (strstr(buf, "SCF Done"))
            break;
    }
    if (sscanf(strchr(buf, '=') + 1, "%lg", & E_n) != 1)
    {
        fprintf(stderr, "Error! Cannot read electron energy of state N! Check your Gaussian output files.");
        Close_file(n_ofl);
        Close_file(np1_ofl);
        Close_file(nm1_ofl);
    }
    /* electron energy of N+1 */
    while (fgets(buf, BUFSIZ, np1_ofl))
    {
        if (strstr(buf, "SCF Done"))
            break;
    }
    if (sscanf(strchr(buf, '=') + 1, "%lg", & E_np1) != 1)
    {
        fprintf(stderr, "Error! Cannot read electron energy of state N+1! Check your Gaussian output files.");
        Close_file(n_ofl);
        Close_file(np1_ofl);
        Close_file(nm1_ofl);
    }
    /* electron energy of N-1 */
    while (fgets(buf, BUFSIZ, nm1_ofl))
    {
        if (strstr(buf, "SCF Done"))
            break;
    }
    if (sscanf(strchr(buf, '=') + 1, "%lg", & E_nm1) != 1)
    {
        fprintf(stderr, "Error! Cannot read electron energy of state N! Check your Gaussian output files.");
        Close_file(n_ofl);
        Close_file(np1_ofl);
        Close_file(nm1_ofl);
    }

    /* HOMO */
    /* there should be at least one line between, so we do not care last_line_pos in the first loop. */
    /* HOMO of N */
    for (;;)
    {
        this_line_pos = ftell(n_ofl);
        if (! fgets(buf, BUFSIZ, n_ofl) || strstr(buf, "Alpha virt."))
            break;
        last_line_pos = this_line_pos;
    }
    fseek(n_ofl, last_line_pos, SEEK_SET);
    fgets(buf, BUFSIZ, n_ofl);
    /* the first string in buf should always be useless. */
    tok = strtok(buf, " \n");
    while ((tok = strtok(NULL, " \n")))
        strcpy(last_value_str, tok);
    if (sscanf(last_value_str, "%lg", & e_HOMO_n) != 1)
    {
        fprintf(stderr, "Error! Cannot read HOMO energy of state N! Check your Gaussian output files.");
        Close_file(n_ofl);
        Close_file(np1_ofl);
        Close_file(nm1_ofl);
    }
    last_line_pos = 0l;
    /* HOMO of N+1 */
    for (;;)
    {
        this_line_pos = ftell(np1_ofl);
        if (! fgets(buf, BUFSIZ, np1_ofl) || strstr(buf, "Alpha virt."))
            break;
        last_line_pos = this_line_pos;
    }
    fseek(np1_ofl, last_line_pos, SEEK_SET);
    fgets(buf, BUFSIZ, np1_ofl);
    tok = strtok(buf, " \n");
    while ((tok = strtok(NULL, " \n")))
        strcpy(last_value_str, tok);
    if (sscanf(last_value_str, "%lg", & e_HOMO_np1) != 1)
    {
        fprintf(stderr, "Error! Cannot read HOMO energy of state N+1! Check your Gaussian output files.");
        Close_file(n_ofl);
        Close_file(np1_ofl);
        Close_file(nm1_ofl);
    }
    
    /* calculates J and J^2 */
    J_n = fabs(e_HOMO_n + E_nm1 - E_n);
    J_np1 = fabs(e_HOMO_np1 + E_n - E_np1);
    * J_ptr = J_n + J_np1;
    * J_squared_ptr = J_n * J_n + J_np1 * J_np1;

    /* close files */
    Close_file(n_ofl);
    Close_file(np1_ofl);
    Close_file(nm1_ofl);

    return;
}

double Calc_J_squared_from_w(double w, void *args)
{
    char sys_command[BUFSIZ + 1] = "";
    FILE *n_ifl = NULL, *np1_ifl = NULL, *nm1_ifl = NULL;
    time_t time_iter_start = 0, time_iter_stop = 0;
    double J = 0.0, J_squared = 0.0;
    char iop_str[] = " IOp(3/107=0000000000,3/108=0000000000)\n"; /* just occupy the position here */
    static unsigned int count_iter = 0u;

    /* prepare files */
    ++ count_iter;
    n_ifl = fopen("N.gjf", "r+t");
    np1_ifl = fopen("Np1.gjf", "r+t");
    nm1_ifl = fopen("Nm1.gjf", "r+t");
    fseek(n_ifl, glob_replace_iop_pos, SEEK_SET);
    fseek(np1_ifl, glob_replace_iop_pos, SEEK_SET);
    fseek(nm1_ifl, glob_replace_iop_pos, SEEK_SET);
    sprintf(iop_str, " IOp(3/107=%05u00000,3/108=%05u00000)\n", \
        (unsigned int)(w * 1E4), (unsigned int)(w * 1E4));
    fprintf(n_ifl, "%s", iop_str);
    fprintf(np1_ifl, "%s", iop_str);
    fprintf(nm1_ifl, "%s", iop_str);
    Close_file(n_ifl);
    Close_file(np1_ifl);
    Close_file(nm1_ifl);
    time_iter_start = time(NULL);

    /* Invoke Gaussian */
    printf("Iteration: %u\n", count_iter);
    printf("w = %6.4lf\n", w);
    sprintf(sys_command, "%s N.gjf N.out", glob_gau_exe);
    printf("Running Gaussian for N state:\n");
    printf("%s\n", sys_command);
    if (system(sys_command))
    {
        fprintf(stderr, "Error! Gaussian job for N state failed.\n");
        fprintf(stderr, "Check your template file and temporary Gaussian output files.\n");
        Print_exit_failure();
    }
    sprintf(sys_command, "%s Np1.gjf Np1.out", glob_gau_exe);
    printf("Running Gaussian for N+1 state:\n");
    printf("%s\n", sys_command);
    if (system(sys_command))
    {
        fprintf(stderr, "Error! Gaussian job for N+1 state failed.\n");
        fprintf(stderr, "Check your template file and temporary Gaussian output files.\n");
        Print_exit_failure();
    }
    sprintf(sys_command, "%s Nm1.gjf Nm1.out", glob_gau_exe);
    printf("Running Gaussian for N-1 state:\n");
    printf("%s\n", sys_command);
    if (system(sys_command))
    {
        fprintf(stderr, "Error! Gaussian job for N-1 state failed.\n");
        fprintf(stderr, "Check your template file and temporary Gaussian output files.\n");
        Print_exit_failure();
    }

    /* Calculates J^2 and J */
    Get_J_and_J_squared(& J, & J_squared);
    time_iter_stop = time(NULL);
    printf("J = %10.8lf, J^2 = %10.8lf\n", J, J_squared);
    printf("Time elapsed for this cycle: %d s.\n", (int)difftime(time_iter_stop, time_iter_start));
    printf("\n");

    return J_squared;
}
