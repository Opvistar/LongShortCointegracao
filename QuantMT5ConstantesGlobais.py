# =============================================================================
#  AUTOR: @Opvistar (Tweeter)
# Arquivo com as constantes globais de configuração 
# =============================================================================

######################################################################
# CONFIGURACAOES LS POR COINTEGRACAO
#######################################################################

#######################################
# CONFIGURACOES DIRETORIOS
PYTHON      = "Database_Python"  # historico de ativos que foram migrados do Tryd
TRYD        = "Database_Tryd"    # historicos de ativos salvos pelo Tryd
MIGRACAO    = "Logs_Migracao"    # logs da migração diária de preços
LONG_SHORT  = "LogsLongShort"    # logs do script long-short
PLANILHAS   = "PlanilhasCSV"     # planilhas geradas para montar as tabelas de controle dos ativos que passaram pelos
                                 # filtros.

#######################################

LISTA_B3 = {'ABEV3','BBAS3','BBDC4','CPLE6','CSNA3','GGBR4','LREN3','PCAR3','TAEE11','VALE3' }

# LISTA_B3 = {'FLRY3','MULT3','MRVE3','CCRO3','QUAL3','LREN3','ABEV3'}
 
STEP_OF_TEST = 10            # steps forward unitroot tests
TAM_AMOSTRA_UNIROOT= 100     # tamanho minimo da amostra para teste uniroot geralmente com n=100

  
AMOSTRA_LONGA = 260          # tamanho MINIMO da amostra para fazer todas os testes, pois o tamanho da amostra do teste uniroot
                             # é para 1 teste, caso queremos correr sobre a série, ela deve ter um tamanho maior.
                             # ATENCAO - SE MUDAR PARA OUTRO VALOR MAIOR QUE 260, VAI EXTRAPOLAR O GRAFICO "*" NR TESTES
                             # MAXIMO É 32.

VALOR_MIN_CORRELA = 0.40     # correlacao minima desejada de coef. correralacao de Pearson, 
                             # ira apresentar os pares com correlacao maior ou igua a VALOR_MIN_CORRELA
MIN_TEST_PASS_RATE= 0.50     # tem que passar em TEST_PASS_RATE dos testes janela movel
                             # diminuiu a volatilidade,,,aumentando pass rate
                         
IS_ADF= False                # adf ou phillips-perron?
OMLY_PASS = False            # só loga os testes que passaram

SPREAD_MINIMO = 1.5

LIMIAR_HORZ = 7.1
HORIZ_FIRST_VALUE = 1
HORIZ_LAST_VALUE  = 0
HORIZ_MEAN_VALUE  = 2 
HORIZ_DISP_RELAT  = 3  
HORIZ_CRITERIO    = HORIZ_DISP_RELAT  

Z_ALFA_0DOT05 = 0.05
Z_ALFA_0DOT01 = 0.01
Z_ALFA_0DOT10 = 0.10

P_VALOR_CRITICO_PP   = Z_ALFA_0DOT05     # se p-valor aumenta, mais dificil rejeitar Ho (não-coestacionariedade)
P_VALOR_CRITICO_KPSS = Z_ALFA_0DOT10     # se p-valor aumenta, rejeita com mais facilidade Ho (estacionariedade)
P_VALOR_CRITICO_REG  = Z_ALFA_0DOT05

 # valor usado para bandas de bollinger, volatilidade etc...geralmente
# é igual ao tamanho do tamanho mínimo dos testes uniroot
N_OTIMO = 160
DISTANCIA_RATIO_BB = 0.01  # razao

# True  = correlacao desabilitada
# False = correlacao minima requerida
COR_OFF = False

