WITH combined_ranks AS (
    SELECT
        individual_id,
        individual_rank_id,
        year,
        group_id,
        rank,
        ord_rank,
        percent_dominated
    FROM individual_ranks
    UNION
    SELECT
        individual_id,
        individual_rank_id,
        year,
        group_id,
        rank,
        ord_rank,
        percent_dominated
    FROM individual_ranks_archive
),
ranks AS (
    SELECT
    ids.individual_code,
    ids.individual_sex,
    ids.individual_id,
    r.year,
    r.group_id,
    g.group_name,
    r.rank,
    r.ord_rank,
    r.percent_dominated
FROM individuals ids
    INNER JOIN combined_ranks r
        ON ids.individual_id = r.individual_id
    INNER JOIN groups g
        ON r.group_id = g.group_id
ORDER BY individual_code, year
)
SELECT
    r.individual_code,
    r.individual_sex,
    r.individual_id,
    r.year,
    r.group_id,
    r.group_name,
    r.rank,
    r.ord_rank,
    r.percent_dominated,
    p.capture_code
FROM ranks r
    INNER JOIN processings p
        ON r.individual_id = p.individual_id;