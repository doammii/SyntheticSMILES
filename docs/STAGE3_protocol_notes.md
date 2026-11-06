# 3단계 — LLaMOL 생성 분자 조건 추가

2단계(2026-10-06 실행)는 `baseline` / `dup` / `randsmi` 세 조건을 끝냈습니다.
3단계는 네 번째 조건 `llamol`을 같은 설계에 추가합니다. 이것으로 비로소
**"생성 분자가 단순 복제를 넘는가"** 에 답할 수 있습니다.

판정 기준은 2단계와 같고, 미리 고정되어 있습니다.

> `llamol`이 `dup`과 `randsmi`를 **둘 다** AUPRC와 MCC에서 이기고,
> seed 5개 중 4개 이상에서 일관되며, assay 2개 이상에서 재현되어야
> 생성모델의 고유 기여로 인정한다. 베이스라인만 이기는 것은 오버샘플링 효과다.

---

## 실행 순서

세 단계이고, **앞 단계 결과를 확인한 뒤** 다음으로 넘어갑니다.

| 순서 | 배치 파일 | 실행 위치 | 예상 시간 |
| --- | --- | --- | --- |
| 1 | `run_stage3_probe.bat` | llamol 저장소 폴더 | 20분 |
| 2 | `run_stage3_generate.bat` | llamol 저장소 폴더 | 3~6시간(추정) |
| 3 | `run_stage3_classify.bat` | `tox21_augmentation` 폴더 | 10~15시간 |

1단계와 2단계는 LLaMOL 저장소 안에서 돌려야 합니다. 스크립트가 저장소의
`sample` 모듈을 import하고 `out\llama2-M-Full-RSS-Canonical.pt`를 읽습니다.
두 배치 파일 맨 위의 `set TOX21=...` 경로를 환경에 맞게 확인하세요.

3단계는 `tox21_augmentation` 폴더에서, `conda activate tox21` 상태로 돌립니다.

---

## 1. 탐침(probe) — 먼저 이것만

`run_stage3_probe.bat`이 두 가지를 확인합니다.

**(1) smoke test** — 조각 `C=CC=O`로 5개를 생성해 모델이 SMILES 컨텍스트를
실제로 반영하는지 봅니다. `fragment retained in 0/5`이 나오면 **거기서 멈추고
알려주세요.** 체크포인트나 저장소 버전 문제이고, 이후 전부가 이 전제에
달려 있습니다.

**(2) 수율 측정** — NR-PPAR-gamma(가장 작은 alert assay)를 실제로 생성해
"원시 생성량 대비 최종 사용 가능 분자 수"를 측정합니다. 이 비율을 모르면
본 생성의 `--per-context` 값이 맞는지 알 수 없습니다.

끝나면 **세 가지를 보고**해 주세요.

- smoke test의 `fragment retained in N/5`
- 마지막 `N unique generated molecules` 숫자
- (2)에 걸린 시간

`N unique generated molecules`가 **372 미만**이면 본 생성 설정이 전부
부족하다는 뜻입니다. 그 경우 알려주시면 값을 다시 계산해 드립니다.

---

## 2. 본 생성

`run_stage3_generate.bat`. assay 5개를 생성합니다.
`--per-context`는 assay마다 다르게 박아두었습니다. alert fragment 수가
적은 assay는 컨텍스트가 적어서 분자당 생성 개수를 크게 올려야 합니다.

| assay | 모드 | fragment | 컨텍스트 | per-context | 필요한 증강 행 |
| --- | --- | --- | --- | --- | --- |
| SR-ARE | alert | 147 (21개 제외) | 1197 | 50 | 1884 |
| SR-MMP | alert | 105 (17개 제외) | 640 | 80 | 1836 |
| SR-p53 | alert | 40 (1개 제외) | 142 | 180 | 846 |
| SR-HSE | alert | 13 (2개 제외) | 129 | 200 | 744 |
| NR-PPAR-gamma | alert | 10 (2개 제외) | 85 | 140 | 372 |
| NR-AR | **whole** | 0 | 308 | 40 | 616 |

**NR-AR은 whole 모드만 가능합니다.** 이 assay는 음성 대조군으로 고른
것이고(BRENK 경고가 활성 분자에서 오히려 고갈, 오즈비 0.70), 그래서
`prepared\fragments\NR-AR.csv`가 설계상 비어 있습니다. alert 조건부 생성이
성립할 수 없으므로 분자 전체를 컨텍스트로 쓰는 arm만 돌립니다.
NR-AR-LBD와 NR-ER-LBD도 같은 이유로 fragment가 0개입니다.

끝나면 `generated\*.csv`의 **행 수**와 콘솔 로그를 보내주세요. 각 파일
마지막에 `WARNING: fewer molecules than the N rows step 2 needs`가 찍혔다면
그 assay는 생성량이 부족한 것이니 3단계로 넘어가기 전에 알려주세요.

---

## 3. 분류 실험

`run_stage3_classify.bat`. 6개 assay × 2개 프로토콜 × `llamol` 조건 ×
seed 5개 = 60회입니다.

결과는 **새 폴더** `out\<assay>_<protocol>_llamol\`에 들어갑니다.
2026-10-06 결과는 그대로 둡니다.

로그에서 **`WARNING:`으로 시작하는 줄을 주의**해 주세요. 그 seed에서
생성 분자가 부족해 추가 행 수가 다른 조건과 맞지 않았다는 뜻이고,
그러면 그 seed는 비교에서 빼야 합니다.

끝나면 `out\*_llamol\results_summary.csv`, `out\*_llamol\results_raw.csv`,
콘솔 로그를 보내주세요. `_work\`(모델 체크포인트)는 용량만 크니 빼주세요.

---

## 2026-10-07 수정 내역

2단계 결과를 받아 3단계 코드를 다시 검토하면서 세 가지 문제를 고쳤습니다.
고치지 않으면 **clean 프로토콜에서 생성 분자가 대부분 또는 전부 버려집니다.**

### (a) 자기 모분자와 매칭되지 않는 fragment

1단계가 뽑은 alert 패턴 일부가 전하 표기 때문에 정작 추출원 분자와
부분구조 매칭이 되지 않습니다(`O=[N+][O-]`, `C=CC(=O)[O]`,
`C=C(C)C(=O)[O]`, `Nc1ccc(N)c([N+])c1` 등).

| assay | 전체 fragment | 매칭 0개 |
| --- | --- | --- |
| SR-ARE | 147 | 21 |
| SR-MMP | 105 | 17 |
| SR-HSE | 13 | 2 |
| NR-PPAR-gamma | 10 | 2 |
| SR-p53 | 40 | 1 |

이전 코드는 이때 `parent`를 fragment 문자열 자체로 대체했습니다
(`parent or row["fragment"]`). 그런 parent는 학습셋에 있을 수 없으므로,
2단계의 "부모가 학습셋에 있는 생성물만 쓴다" 필터를 절대 통과하지 못하고
해당 fragment에서 나온 생성물이 **조용히 전량 폐기**됩니다.
지금은 이런 fragment를 건너뛰고 개수를 로그에 찍습니다.

### (b) fragment당 parent가 1개

이전에는 fragment를 포함하는 **첫 번째** 활성 분자만 parent로 썼습니다.
그 분자가 해당 seed의 valid/test에 들어가면 그 fragment의 생성물 블록이
통째로 사라집니다. SR-HSE는 고유 parent가 39개, NR-PPAR-gamma는 30개뿐이라
seed마다 사용 가능 수량이 크게 흔들립니다.

지금은 `--parents-per-fragment`(기본 25)만큼, 그 fragment를 실제로 포함하는
서로 다른 활성 분자를 섞어서 짝지웁니다. parent가 분할 전체에 퍼지고,
fragment가 적은 assay의 컨텍스트 수도 함께 늘어납니다.

실측하면 이전 코드의 상태는 추정보다 나빴습니다. 아래 "이전"은 fragment당
parent 1개 배정에서 **실제 분자였던 고유 parent 수**이고, 나머지는 (a)의
fragment 문자열이라 생성물이 전량 폐기됩니다.

| assay | 이전 (실제 분자) | 수정 후 | 컨텍스트 |
| --- | --- | --- | --- |
| SR-ARE | 73 | 389 | 1197 |
| SR-MMP | 50 | 271 | 640 |
| SR-p53 | 21 | 44 | 142 |
| SR-HSE | **5** | 39 | 129 |
| NR-PPAR-gamma | **4** | 30 | 85 |

SR-HSE는 생성 분자 744개를 **5개 분자**에서 파생시키는 셈이었습니다.
화학공간 분석을 하면 다양성이 바닥으로 나올 수밖에 없는 구성입니다.

검증: 5개 seed에서 "부모가 학습셋에 있는 생성물" 비율이 평균 63~66%로
안정적이고, **테스트셋에 들어간 생성 분자는 모든 seed에서 0개**입니다.

### (c) `--max-contexts`의 조용한 절단

기본값 400으로 컨텍스트 목록을 말없이 잘랐습니다. SR-ARE(1197)와
SR-MMP(640)가 영향을 받습니다. 기본값을 0(무제한)으로 바꿨습니다.

### (d) 2단계 `--n-add` — 배치 파일에 고정

`02_run_experiments.py`는 추가 행 수를 `--synthetic`이 없으면
`2 × 소수클래스`로, **있으면 `len(synthetic)`으로** 정합니다.
그대로 돌리면 추가 행 수가 바뀌어 2026-10-06 결과와 비교할 수 없습니다.
`run_stage3_classify.bat`에 assay별 값(1884 / 1836 / 846 / 744 / 372 / 616)을
명시해 두었으니, 직접 명령을 입력할 때도 `--n-add`를 빼지 마세요.

`02_run_experiments.py`에는 생성 분자가 부족해 추가 행 수가 맞지 않을 때
경고를 찍는 줄만 추가했습니다. 계산 로직은 건드리지 않았으므로
2026-10-06 결과는 그대로 재현됩니다.
